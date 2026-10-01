import json
import re
import socket
import urllib.request
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Optional
from urllib.error import HTTPError, URLError

import numpy as np
from fastapi import HTTPException

from app.schemas.topology import TMParams, TopologyRegion
from app.services.topology.labels import (
    INTRA, SIGNAL, TM, UNASSIGNED, region_fields, side_from_text,
)
from app.services.topology.providers.tm.base import TMPrediction, TMProvider
from app.services.topology.residues import ResidueFrame, load_residue_frame, map_to_reference

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{acc}.json"
# Official UniProtKB accession format (+ optional isoform suffix). Validating it also
# keeps user input from being spliced into the request path.
ACCESSION_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?$")
FEATURE_TYPES = ("Topological domain", "Signal", "Intramembrane", "Transmembrane")
PAINT_ORDER = {t: i for i, t in enumerate(FEATURE_TYPES)}   # later types win on overlap


@dataclass
class UniProtReference:
    accession: str
    offset: Optional[int] = None      # UniProt position - author resseq, when the file says so


def _as_list(value):
    if value is None:
        return []
    return list(value) if isinstance(value, (list, tuple)) else [value]


def _int_or_none(text):
    try:
        return int(str(text).strip())
    except (TypeError, ValueError):
        return None


def _reference_from_pdb(file_path: Path, chain_id: Optional[str]) -> Optional[UniProtReference]:
    """DBREF / DBREF1+DBREF2 records (Bio.PDB's header parser does not read them)."""
    found = []
    pending = {}
    with open(file_path, errors="replace") as fh:
        for line in fh:
            rec = line[:6]
            if rec == "DBREF ":
                if line[26:32].strip() != "UNP":
                    continue
                found.append((line[12:13].strip(), line[33:41].strip(),
                              _int_or_none(line[14:18]), _int_or_none(line[55:60])))
            elif rec == "DBREF1" and line[26:32].strip() == "UNP":
                pending[line[12:13].strip()] = _int_or_none(line[14:18])
            elif rec == "DBREF2":
                ch = line[12:13].strip()
                if ch in pending:
                    found.append((ch, line[18:40].strip(), pending.pop(ch), _int_or_none(line[45:55])))
            elif rec in ("ATOM  ", "HETATM"):
                break                       # header is over
    if not found:
        return None
    chosen = next((f for f in found if f[0] == chain_id), found[0])
    _, acc, seq_begin, db_begin = chosen
    offset = db_begin - seq_begin if seq_begin is not None and db_begin is not None else None
    return UniProtReference(acc, offset)


def _reference_from_cif(file_path: Path, chain_id: Optional[str]) -> Optional[UniProtReference]:
    """_struct_ref (+ _struct_ref_seq for the chain / offset); AlphaFold's
    _ma_target_ref_db_details as fallback. (MMCIFParser has no get_dict(): use MMCIF2Dict.)"""
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    d = MMCIF2Dict(str(file_path))
    ref_ids = _as_list(d.get("_struct_ref.id"))
    db_names = _as_list(d.get("_struct_ref.db_name"))
    accessions = _as_list(d.get("_struct_ref.pdbx_db_accession"))
    unp = {rid: acc for rid, db, acc in zip(ref_ids, db_names, accessions) if db == "UNP"}
    if unp:
        seq_ref = _as_list(d.get("_struct_ref_seq.ref_id"))
        strands = _as_list(d.get("_struct_ref_seq.pdbx_strand_id"))
        auth_beg = _as_list(d.get("_struct_ref_seq.pdbx_auth_seq_align_beg"))
        db_beg = _as_list(d.get("_struct_ref_seq.db_align_beg"))
        rows = list(zip(seq_ref, strands, auth_beg, db_beg))
        for rid, strand, a_beg, u_beg in rows:
            if rid in unp and (chain_id is None or strand == chain_id):
                a, u = _int_or_none(a_beg), _int_or_none(u_beg)
                return UniProtReference(unp[rid], u - a if a is not None and u is not None else None)
        return UniProtReference(next(iter(unp.values())))
    af_db = _as_list(d.get("_ma_target_ref_db_details.db_name"))
    af_acc = _as_list(d.get("_ma_target_ref_db_details.db_accession"))
    for db, acc in zip(af_db, af_acc):
        if db == "UNP":
            return UniProtReference(acc, 0)
    return None


def extract_uniprot_reference(file_path: Path, chain_id: Optional[str] = None) -> Optional[UniProtReference]:
    try:
        if file_path.suffix.lower() in (".cif", ".mmcif"):
            return _reference_from_cif(file_path, chain_id)
        return _reference_from_pdb(file_path, chain_id)
    except (OSError, ValueError, KeyError) as err:
        print(f"[uniprot] could not read UniProt reference from {file_path.name}: {err}")
        return None


@lru_cache(maxsize=256)
def fetch_uniprot_entry(accession: str) -> dict:
    """UniProt REST JSON. Cached per accession (errors are not cached)."""
    req = urllib.request.Request(UNIPROT_URL.format(acc=accession),
                                 headers={"User-Agent": "ProteinVisualizeApp/1.0",
                                          "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=10) as response:
            payload = response.read()
    except HTTPError as err:
        if err.code in (400, 404):
            raise HTTPException(status_code=404, detail=f"UniProt ID {accession} not found") from err
        raise HTTPException(status_code=502, detail=f"UniProt API error: {err.reason}") from err
    except (socket.timeout, TimeoutError) as err:        # read timeouts are not URLError
        raise HTTPException(status_code=504, detail="UniProt API timed out") from err
    except URLError as err:
        raise HTTPException(status_code=504,
                            detail=f"Network error connecting to UniProt API: {err.reason}") from err
    try:
        return json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as err:
        raise HTTPException(status_code=502, detail="UniProt API returned invalid JSON") from err


_TM_NAME_RE = re.compile(r"Name=([^;\s]+)")
MIN_REPORTED_RUN = 3        # numbering runs shorter than this (insertion codes...) are not reported


def _tm_name(feature: dict, index: int) -> str:
    """'Helical; Name=17' -> 'TM17'; unnamed helices get their order in the entry."""
    m = _TM_NAME_RE.search(feature.get("description") or "")
    return f"TM{m.group(1)}" if m else f"TM{index}"


def _compress_names(names: list[str]) -> str:
    """['TM1', 'TM2', 'TM3', 'TM7'] -> 'TM1-TM3, TM7' (only numeric names are joined)."""
    out: list[str] = []
    run: list[int] = []

    def flush():
        if run:
            out.append(f"TM{run[0]}" if len(run) == 1 else f"TM{run[0]}-TM{run[-1]}")
            run.clear()

    for name in names:
        num = name[2:]
        if num.isdigit() and run and int(num) == run[-1] + 1:
            run.append(int(num))
        elif num.isdigit():
            flush()
            run.append(int(num))
        else:
            flush()
            out.append(name)
    flush()
    return ", ".join(out)


def _numbering_runs(frame: ResidueFrame, uniprot_pos: list[Optional[int]]):
    """Stretches where (UniProt position - author resseq) is constant:
    [(delta, first frame position, last frame position)]."""
    runs: list[list[int]] = []
    for p, u in enumerate(uniprot_pos):
        if u is None:
            continue
        d = u - frame.residues[p].resseq
        if runs and runs[-1][0] == d:
            runs[-1][2] = p
        else:
            runs.append([d, p, p])
    return [tuple(r) for r in runs]


def _feature_bounds(feature: dict) -> Optional[tuple[int, int]]:
    loc = feature.get("location") or {}
    start = _int_or_none((loc.get("start") or {}).get("value"))
    end = _int_or_none((loc.get("end") or {}).get("value"))
    if start is None or end is None:
        return None                 # unknown / uncertain boundary ("?")
    return (start, end) if start <= end else (end, start)


class UniprotTMProvider(TMProvider):
    def _extract_uniprot_id_from_pdb(self, file_path: Path, chain_id: Optional[str] = None) -> str:
        """Kept for older callers."""
        ref = extract_uniprot_reference(Path(file_path), chain_id)
        return ref.accession if ref else ""

    def predict_tm(self, file_path: Path, params: Optional[TMParams] = None,
                   frame: Optional[ResidueFrame] = None, **kwargs) -> TMPrediction:
        frame = frame if frame is not None else load_residue_frame(file_path, kwargs.get("chain_id"))
        warns: list[str] = []

        ref = None
        given = (kwargs.get("uniprot_id") or "").strip().upper()
        if given:
            ref = UniProtReference(given)
            file_ref = extract_uniprot_reference(file_path, frame.chain_id)
            if file_ref and file_ref.accession.upper() == given:
                ref.offset = file_ref.offset
        else:
            ref = extract_uniprot_reference(file_path, frame.chain_id)
        if not ref or not ref.accession:
            raise HTTPException(
                status_code=400,
                detail="UniProt ID is required for the UniProt algorithm. Could not extract it "
                       "from the structure file. Please provide it manually.")
        acc = ref.accession.strip().upper()
        if not ACCESSION_RE.match(acc):
            raise HTTPException(status_code=400, detail=f"'{acc}' is not a valid UniProt accession")

        data = fetch_uniprot_entry(acc)
        useq = (data.get("sequence") or {}).get("value", "")
        if not useq:
            raise HTTPException(status_code=502, detail=f"UniProt entry {acc} has no sequence")
        if len(frame) == 0:
            return TMPrediction(labeler=f"UniProt_{acc}", warnings=["no amino-acid residues in chain"])

        # UniProt position u (1-based) <-> frame position, verified by residue identity
        offsets = [0] + ([ref.offset] if ref.offset else [])
        mapping = map_to_reference(frame, [(u, "") for u in range(1, len(useq) + 1)],
                                   list(useq), offsets=offsets, source="UniProt",
                                   allow_order_fallback=False, note_remap=False)
        if mapping.matched == 0:
            raise HTTPException(
                status_code=422,
                detail=f"UniProt {acc} sequence does not match chain {frame.chain_id} of the "
                       f"structure (identity {mapping.identity or 0:.0%}). Wrong UniProt ID or chain?")
        warns.extend(mapping.warnings)
        if mapping.method != "residue numbering":
            warns.append(f"UniProt {acc}: UniProt positions mapped to chain {frame.chain_id} by "
                         f"{mapping.method} (identity {mapping.identity or 0:.0%}, "
                         f"{mapping.matched}/{len(frame)} residues)")

        pos_of_u = np.full(len(useq) + 2, -1, dtype=int)
        uniprot_pos: list[Optional[int]] = [None] * len(frame)
        for p, j in enumerate(mapping.ref_index):
            if j is not None:
                pos_of_u[j + 1] = p
                uniprot_pos[p] = j + 1

        # Author numbering is what every region of the result is reported in. When it does
        # not follow UniProt (8ZU3 / 9VMX: author 767-892 = UniProt 789-914), say so - a
        # TM shown at 795-806 is UniProt TM17 (817-828), not a mapping error.
        runs = _numbering_runs(frame, uniprot_pos)
        shifted = [(d, a, b) for d, a, b in runs if d != 0 and b - a + 1 >= MIN_REPORTED_RUN]
        if shifted:
            parts = [f"author {frame.residues[a].label}-{frame.residues[b].label} = UniProt "
                     f"{uniprot_pos[a]}-{uniprot_pos[b]} ({d:+d})" for d, a, b in shifted[:6]]
            more = f" and {len(shifted) - 6} more stretch(es)" if len(shifted) > 6 else ""
            warns.append(f"UniProt {acc}: chain {frame.chain_id} author residue numbers differ "
                         f"from UniProt numbering - {'; '.join(parts)}{more}. Residues are "
                         "matched by sequence alignment; positions in this result are AUTHOR "
                         "numbers (the TM names below give the UniProt range)")

        n = len(frame)
        labels = [UNASSIGNED] * n
        segments: list[tuple[int, int]] = []
        display: list[tuple[int, TopologyRegion]] = []
        outside = 0
        tm_total = 0
        tm_missing: list[str] = []                    # not resolved at all -> skipped
        tm_partial: list[str] = []                    # only partly resolved -> kept
        tm_shifted: list[str] = []                    # resolved, author numbers != UniProt
        features = [f for f in data.get("features", []) if f.get("type") in FEATURE_TYPES]
        tm_names = {id(f): _tm_name(f, k + 1) for k, f in enumerate(sorted(
            (f for f in features if f["type"] == "Transmembrane"),
            key=lambda f: (_feature_bounds(f) or (0, 0))[0]))}
        features.sort(key=lambda f: (PAINT_ORDER[f["type"]], (_feature_bounds(f) or (0, 0))[0]))
        for f in features:
            bounds = _feature_bounds(f)
            if bounds is None:
                continue
            ftype = f["type"]
            u0, u1 = max(1, bounds[0]), min(len(useq), bounds[1])
            hits = pos_of_u[u0:u1 + 1] if u0 <= u1 else np.array([], dtype=int)
            hits = hits[hits >= 0]
            if ftype == "Transmembrane":
                tm_total += 1
            if hits.size == 0:
                outside += 1          # feature not resolved in this model: do NOT snap it
                if ftype == "Transmembrane":
                    tm_missing.append(tm_names[id(f)])
                continue
            lo, hi = int(hits.min()), int(hits.max())
            desc = f.get("description") or ftype
            if ftype == "Transmembrane":
                name = tm_names[id(f)]
                desc = (f"{desc} (UniProt {u0}-{u1})" if "Name=" in desc
                        else f"{desc} ({name}, UniProt {u0}-{u1})")
                label = TM
                segments.append((lo, hi))
                if hits.size < u1 - u0 + 1:
                    tm_partial.append(f"{name} {u0}-{u1}: {hits.size}/{u1 - u0 + 1} resolved")
                a_lo, a_hi = frame.residues[lo], frame.residues[hi]
                if (a_lo.resseq, a_hi.resseq) != (uniprot_pos[lo], uniprot_pos[hi]):
                    tm_shifted.append(f"{name} {uniprot_pos[lo]}-{uniprot_pos[hi]} = author "
                                      f"{a_lo.label}-{a_hi.label}")
            elif ftype == "Intramembrane":
                label = INTRA
            elif ftype == "Signal":
                label = SIGNAL
            else:
                label = side_from_text(desc) or UNASSIGNED
                if label == UNASSIGNED:
                    warns.append(f"UniProt topological domain '{desc}' has no known side; "
                                 "side inferred from topology")
            for k in range(lo, hi + 1):
                labels[k] = label
            _, side, _ = region_fields(label)
            display.append((lo, TopologyRegion(
                type=ftype, start=frame.residues[lo].resseq, end=frame.residues[hi].resseq,
                description=desc, side=side if ftype != "Transmembrane" else "membrane",
                start_icode=frame.residues[lo].icode or None,
                end_icode=frame.residues[hi].icode or None)))
        if tm_total:
            msg = (f"UniProt {acc}: {tm_total - len(tm_missing)} of {tm_total} UniProt TM helices "
                   f"are resolved in chain {frame.chain_id}")
            if tm_missing:
                msg += (f"; not resolved (no coordinates, skipped): "
                        f"{_compress_names(tm_missing)}")
            warns.append(msg)
        other_missing = outside - len(tm_missing)
        if other_missing:
            warns.append(f"UniProt {acc}: {other_missing} other topology feature(s) lie outside "
                         "the residues resolved in this structure and were skipped")
        if tm_partial:
            warns.append(f"UniProt {acc}: partly resolved TM helices, drawn over their resolved "
                         f"residues only - {'; '.join(tm_partial)}")
        if tm_shifted:
            warns.append(f"UniProt {acc}: TM helices whose author numbers differ from UniProt - "
                         f"{'; '.join(tm_shifted)}")

        display.sort(key=lambda item: item[0])
        return TMPrediction(
            regions=[r for _, r in display],
            labels=labels,
            segments=sorted(segments),
            labeler=f"UniProt_{acc}",
            warnings=warns,
        )
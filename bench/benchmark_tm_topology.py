"""
Reproducible benchmark of the ProteinLab transmembrane-topology engine.

Ground truth
------------
UniProtKB/Swiss-Prot "Transmembrane" features.  UniProt positions are mapped onto
the analysed chain of the PDB entry through the *same* residue-alignment machinery
the application uses at run time (`app.services.topology.residues.map_to_reference`),
so a difference in a reported boundary is a difference in the *prediction*, not in the
numbering.  A reference helix enters the benchmark only when at least
``MIN_RESOLVED_FRACTION`` of its residues have coordinates in the analysed chain.

Predicted methods
-----------------
3d_energy        implicit-membrane transfer-energy placement + bilayer labelling
3d_slab_geom     500-axis signed-hydropathy slab fit + labelling
kyte_doolittle_seq  1-D two-scale hydropathy scan
tmhmm            TMHMM 2.0 (vendored numpy port of pyTMHMM)
uniprot_api      the curated UniProt annotation itself, painted on the structure

Metrics
-------
One-to-one greedy matching of predicted segments to reference segments by overlap.
A pair is *matched* when the overlap covers at least ``MATCH_IOU`` of the smaller
segment.  For matched pairs the two boundary errors (|d start|, |d end|) are
reported in residues.

Usage
-----
    python bench/benchmark_tm_topology.py --out bench/results --limit 60
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import random
import re
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Optional

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))
sys.path.insert(0, str(ROOT / "bench"))

import numpy as np  # noqa: E402

from app.schemas.topology import TMParams  # noqa: E402
from app.services.topology.providers.tm.base import TMPrediction  # noqa: E402
from app.services.topology.providers.tm.energy import EnergyTMProvider  # noqa: E402
from app.services.topology.providers.tm.sequence import SequenceTMProvider  # noqa: E402
from app.services.topology.providers.tm.sequence import SequenceTMProvider  # noqa: E402
from app.services.topology.providers.tm.tmhmm import TMHMMProvider  # noqa: E402
from app.services.topology.providers.tm.uniprot import (  # noqa: E402
    UniprotTMProvider, extract_uniprot_reference, fetch_uniprot_entry,
)
from app.services.topology.residues import load_residue_frame, map_to_reference  # noqa: E402
from app.services.topology.topology_predictor import GeometryTMProvider  # noqa: E402

UA = {"User-Agent": "ProteinLab-benchmark/1.0 (academic validation)"}

# ---- benchmark configuration (fixed -> the benchmark is reproducible) ----------
SEED = 20261021
MIN_TM_ANNOTATIONS = 3          # proteins with fewer curated helices are too small
MAX_TM_ANNOTATIONS = 16
MIN_RESOLVED_FRACTION = 0.60    # a curated helix must be >=60% resolved in the chain
MATCH_IOU = 0.50                # overlap fraction of the smaller segment
BOUNDARY_TOLERANCE = 5          # residues, used for the "within tolerance" summary
SOLUBLE_QUERY = "(reviewed:true) AND (ft_transmem:*) AND (organism_id:9606)"

# A fixed, cited control set of soluble proteins (Swiss-Prot accessions); the PDB
# entry and chain are chosen by the script (best resolution, first chain).
SOLUBLE_ACCESSIONS = [
    "P00253",   # Lysozyme C, chicken
    "P00698",   # Lysozyme C precursor, chicken
    "P0CG47",   # Polyubiquitin-B, human
    "P02185",   # Myoglobin, human
    "P00099",   # Cytochrome c, human
    "P42212",   # Green fluorescent protein, Aequorea victoria
    "P02752",   # Fibronectin, human
    "P02766",   # Transthyretin, human
    "P04040",   # Catalase, human
    "P69441",   # Adenylate kinase, E. coli
    "P01568",   # Crambin, cress
    "P00760",   # Trypsin, bovine
    "P0DP23",   # Calmodulin-1, human
    "P00533",   # Epidermal growth factor receptor kinase domain, human (soluble
                #   fragment of a membrane protein - a deliberately hard control)
    "P10809",   # 60 kDa heat shock protein, human (large soluble oligomer)
]

METHODS = {
    "3d_energy": EnergyTMProvider,
    "3d_slab_geom": GeometryTMProvider,
    "kyte_doolittle_seq": SequenceTMProvider,
    "tmhmm": TMHMMProvider,
    "uniprot_api": UniprotTMProvider,
}


# --------------------------------------------------------------------------- io
def http_json(url: str, timeout: int = 60):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def download(url: str, dest: Path, timeout: int = 120) -> bool:
    if dest.exists() and dest.stat().st_size > 1000:
        return True
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = resp.read()
    except Exception:
        return False
    if len(data) < 1000:
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(data)
    return True


def uniprot_search(query: str, size: int = 500) -> tuple[list[tuple[str, list[str]]], Optional[str]]:
    """Return [(accession, [pdb_id, ...])] for the query.

    UniProt's search language has no valid cross-reference filter for PDB, so the
    `xref_pdb` column is requested and filtered here.
    """
    out: list[tuple[str, list[str]]] = []
    cursor = None
    release = None
    while len(out) < size:
        params = {"query": query, "format": "tsv", "fields": "accession,xref_pdb",
                  "size": "500"}
        if cursor:
            params["cursor"] = cursor
        url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(params)
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=90) as resp:
            body = resp.read().decode("utf-8")
            # the cursor of the next page is exposed through the RFC 5988 Link header
            link = resp.headers.get("Link") or ""
            release = resp.headers.get("X-UniProt-Release") or release
            next_cursor = None
            if 'rel="next"' in link:
                m = re.search(r"cursor=([A-Za-z0-9]+)", link)
                if m:
                    next_cursor = m.group(1)
        lines = body.strip().splitlines()
        for ln in lines[1:]:
            parts = ln.split("\t")
            if len(parts) < 2:
                continue
            acc = parts[0].strip()
            pdbs = [p.strip() for p in parts[1].split(",") if p.strip()]
            if acc and pdbs:
                out.append((acc, pdbs))
        if not next_cursor:
            break
        cursor = next_cursor
    return out[:size], release


# ------------------------------------------------------------------ ground truth
def reference_segments(frame, entry: dict, accession: str):
    """UniProt 'Transmembrane' features as (start_pos, end_pos) in the frame."""
    useq = (entry.get("sequence") or {}).get("value", "")
    if not useq:
        return [], 0
    ref = extract_uniprot_reference(Path(frame.source_path), frame.chain_id)
    offset = ref.offset if ref and ref.accession.upper() == accession.upper() else None
    offsets = [0] + ([offset] if offset else [])
    mapping = map_to_reference(frame, [(u, "") for u in range(1, len(useq) + 1)],
                              list(useq), offsets=offsets, source="UniProt",
                              allow_order_fallback=False, note_remap=False)
    if mapping.matched == 0:
        return [], 0
    pos_of_u = np.full(len(useq) + 2, -1, dtype=int)
    for p, j in enumerate(mapping.ref_index):
        if j is not None:
            pos_of_u[j + 1] = p

    out = []
    for feat in entry.get("features", []):
        if feat.get("type") != "Transmembrane":
            continue
        loc = feat.get("location") or {}
        u0 = (loc.get("start") or {}).get("value")
        u1 = (loc.get("end") or {}).get("value")
        try:
            u0, u1 = int(u0), int(u1)
        except (TypeError, ValueError):
            continue
        if u0 > u1:
            u0, u1 = u1, u0
        u0, u1 = max(1, u0), min(len(useq), u1)
        hits = pos_of_u[u0:u1 + 1]
        hits = hits[hits >= 0]
        if hits.size == 0:
            continue
        if hits.size < MIN_RESOLVED_FRACTION * (u1 - u0 + 1):
            continue
        out.append((int(hits.min()), int(hits.max())))
    return sorted(out), int(mapping.matched)


# ---------------------------------------------------------------------- metrics
def overlap(a, b) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


def match_segments(ref, pred):
    """Greedy one-to-one matching by overlap. Returns (pairs, n_matched, mean_dstart...)."""
    pairs = []
    used = set()
    for i, r in enumerate(ref):
        best_j, best_ov = None, 0
        for j, p in enumerate(pred):
            if j in used:
                continue
            ov = overlap(r, p)
            if ov > best_ov:
                best_j, best_ov = j, ov
        if best_j is None:
            continue
        r_len = r[1] - r[0] + 1
        p_len = pred[best_j][1] - pred[best_j][0] + 1
        if best_ov >= MATCH_IOU * min(r_len, p_len):
            used.add(best_j)
            pairs.append((i, best_j, best_ov))
    return pairs, len(pairs)


def score(ref, pred):
    pairs, n_match = match_segments(ref, pred)
    ds, de = [], []
    for i, j, _ in pairs:
        ds.append(abs(pred[j][0] - ref[i][0]))
        de.append(abs(pred[j][1] - ref[i][1]))
    ref_res = set()
    for s, e in ref:
        ref_res.update(range(s, e + 1))
    hit = 0
    for s, e in pred:
        hit += len(set(range(s, e + 1)) & ref_res)
    recall = n_match / len(ref) if ref else 0.0
    precision = n_match / len(pred) if pred else 0.0
    f1 = 2 * recall * precision / (recall + precision) if (recall + precision) else 0.0
    both = [max(a, b) for a, b in zip(ds, de)]
    return {
        "n_ref": len(ref), "n_pred": len(pred), "n_match": n_match,
        "recall": round(recall, 4), "precision": round(precision, 4), "f1": round(f1, 4),
        "residue_recall": round(hit / len(ref_res), 4) if ref_res else 0.0,
        "mean_dstart": round(float(np.mean(ds)), 2) if ds else None,
        "mean_dend": round(float(np.mean(de)), 2) if de else None,
        "median_max_boundary": float(np.median(both)) if both else None,
        "within_tol_frac": round(sum(1 for x in both if x <= BOUNDARY_TOLERANCE) / len(both), 4)
        if both else None,
    }


# ------------------------------------------------------------------------- main
def pick_pdb(entry: dict):
    """Return (pdb_id, chain) for a UniProt entry from its PDB cross-references.

    The entry with the *best* (lowest) resolution is preferred.  UniProt writes the
    resolution as e.g. "1.70 A" and "-" for NMR and predicted models, and names the
    chains as e.g. "A=2-149" or "B/D/F/H=1-149"; both are parsed here.
    """
    best = None
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") != "PDB":
            continue
        props = {p.get("key"): p.get("value") for p in xref.get("properties", [])}
        pdb_id = xref.get("id")
        chain = _first_chain(props.get("Chains") or props.get("Chain") or "")
        res = _resolution(props.get("Resolution"))
        if best is None or res < best[0]:
            best = (res, pdb_id, chain)
    if best is None:
        return None, None
    return best[1], best[2]


def _resolution(raw) -> float:
    """'1.70 A' -> 1.70; '-' or missing -> 99.0 (ranked last)."""
    if raw is None:
        return 99.0
    m = re.match(r"\s*([0-9]+(?:\.[0-9]+)?)", str(raw))
    return float(m.group(1)) if m else 99.0


def _first_chain(raw) -> str:
    """'B/D/F/H=1-149' -> 'B';  'A=2-149' -> 'A';  '' -> ''."""
    if not raw:
        return ""
    return re.split(r"[,/=]", str(raw).strip())[0].strip()


def analyse_one(job) -> Optional[dict]:
    """Download, parse and score one UniProt accession. Returns None when the entry
    does not qualify (no curated helices, no PDB cross-reference, too few resolved
    helices).  Module level so that it can be pickled by ProcessPoolExecutor."""
    acc, structs_dir = job
    structs = Path(structs_dir)
    try:
        entry = fetch_uniprot_entry(acc)
    except Exception as exc:
        print(f"[bench] {acc}: UniProt fetch failed ({exc})")
        return None
    tm_feats = [f for f in entry.get("features", []) if f.get("type") == "Transmembrane"]
    if not (MIN_TM_ANNOTATIONS <= len(tm_feats) <= MAX_TM_ANNOTATIONS):
        return None
    pdb_id, chain = pick_pdb(entry)
    if not pdb_id:
        return None
    path = structs / f"{pdb_id}.cif"
    if not download(f"https://files.rcsb.org/download/{pdb_id.upper()}.cif", path):
        print(f"[bench] {acc} / {pdb_id}: download failed")
        return None

    try:
        frame = load_residue_frame(path, chain or None)
    except Exception as exc:
        print(f"[bench] {acc} / {pdb_id}: parse failed ({exc})")
        return None
    if len(frame) < 50:
        return None
    # the analysed chain must be the one that maps onto the UniProt sequence
    ref, matched = reference_segments(frame, entry, acc)
    if len(ref) < 2:
        return None

    row = {"accession": acc, "pdb_id": pdb_id, "chain": frame.chain_id,
           "protein_name": (entry.get("proteinDescription") or {})
           .get("recommendedName", {}).get("fullName", {}).get("value", ""),
           "keywords": entry.get("keywords", []),
           "n_residues": len(frame), "n_reference_tm": len(ref),
           "mapping_matched": matched, "methods": {}}
    for name, cls in METHODS.items():
        t0 = time.time()
        try:
            pred: TMPrediction = cls().predict_tm(path, params=TMParams(), frame=frame)
            segs = [(int(a), int(b)) for a, b in (pred.segments or [])]
            sc = score(ref, segs)
            sc["seconds"] = round(time.time() - t0, 2)
            if name == "3d_energy":
                sc["dG_transfer"] = round(-float(pred.membrane_score), 2)
            row["methods"][name] = sc
        except Exception as exc:
            row["methods"][name] = {"error": f"{type(exc).__name__}: {exc}"}

    # --- reference predictors: published rules implemented for the comparison only
    t_refs = time.time()
    try:
        import reference_predictors as RP

        kd = [(int(a), int(b)) for a, b in
              (SequenceTMProvider().predict_tm(path, params=TMParams(), frame=frame).segments or [])]
        refs = RP.reference_predictors(frame.sequence, path, len(frame),
                                       hydropathy_segments=kd, chain_id=frame.chain_id)
        for name, segs in refs.items():
            sc = score(ref, [(int(a), int(b)) for a, b in segs])
            sc["seconds"] = round(time.time() - t_refs, 2)
            row["methods"][name] = sc
    except Exception as exc:
        row.setdefault("reference_error", f"{type(exc).__name__}: {exc}")
    return row


def run(out_dir: Path, limit: int, query: str, seed: int = SEED, workers: int = 1):
    out_dir.mkdir(parents=True, exist_ok=True)
    structs = out_dir / "structures"
    results = []
    meta = {"query": query, "seed": seed, "limit": limit, "workers": workers,
            "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "uniprot_release": None, "match_iou": MATCH_IOU,
            "min_resolved_fraction": MIN_RESOLVED_FRACTION,
            "boundary_tolerance_residues": BOUNDARY_TOLERANCE}

    accessions, release = uniprot_search(query, size=4000)
    random.Random(seed).shuffle(accessions)
    meta["uniprot_release"] = release
    meta["n_candidates"] = len(accessions)
    print(f"[bench] {len(accessions)} candidate accessions with curated TM helices "
          f"and a PDB cross-reference (UniProt release {release})")

    jobs = [(acc, str(structs)) for acc, _pdbs in accessions]

    # the cohort is a deterministic prefix of the shuffled candidate list; entries that
    # do not qualify are skipped, so the scan continues until `limit` proteins qualify
    n_done = 0
    if workers > 1:
        pool = ProcessPoolExecutor(max_workers=workers)
        try:
            it = iter(jobs)
            window = workers * 2
            pending = []
            for job in itertools.islice(it, window):
                pending.append(pool.submit(analyse_one, job))
            exhausted = False
            while pending and n_done < limit:
                future = pending.pop(0)
                row = future.result()
                if row is not None:
                    n_done += 1
                    results.append(row)
                    print(f"[bench] {n_done:3d} {row['accession']} {row['pdb_id']} "
                          f"{row['chain']} refTM={row['n_reference_tm']} "
                          + " ".join(f"{k.split('_')[0]}={row['methods'][k].get('f1', '-')}"
                                     for k in row["methods"] if "error" not in row["methods"][k]))
                    (out_dir / "benchmark.json").write_text(json.dumps(results, indent=2))
                if not exhausted:
                    try:
                        pending.append(pool.submit(analyse_one, next(it)))
                    except StopIteration:
                        exhausted = True
        finally:
            # wait=True and cancel_futures keep the worker processes from being
            # orphaned on Windows when the script is interrupted
            pool.shutdown(wait=True, cancel_futures=True)
    else:
        for job in jobs:
            row = analyse_one(job)
            if row is None:
                continue
            n_done += 1
            results.append(row)
            print(f"[bench] {n_done:3d} {row['accession']} {row['pdb_id']} {row['chain']} "
                  f"refTM={row['n_reference_tm']} "
                  + " ".join(f"{k.split('_')[0]}={row['methods'][k].get('f1', '-')}"
                             for k in row["methods"] if "error" not in row["methods"][k]))
            (out_dir / "benchmark.json").write_text(json.dumps(results, indent=2))
            if n_done >= limit:
                break
    meta["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    meta["n_proteins"] = len(results)
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2))
    return results


def run_soluble(out_dir: Path):
    """Control cohort: the transfer free energy of a soluble protein must not look
    like a membrane.  A structure is called 'membrane' when dG_transfer <= -10 kcal/mol."""
    out_dir.mkdir(parents=True, exist_ok=True)
    structs = out_dir / "structures"
    rows = []
    for acc in SOLUBLE_ACCESSIONS:
        try:
            entry = fetch_uniprot_entry(acc)
        except Exception as exc:
            print(f"[bench] soluble {acc}: UniProt fetch failed ({exc})")
            continue
        pdb_id, chain = pick_pdb(entry)
        if not pdb_id:
            continue
        path = structs / f"{pdb_id}.cif"
        if not download(f"https://files.rcsb.org/download/{pdb_id.upper()}.cif", path):
            continue
        try:
            frame = load_residue_frame(path, chain or None)
        except Exception:
            continue
        if len(frame) < 30:
            continue
        dg, thickness, curved = None, None, None
        try:
            pred = EnergyTMProvider().predict_tm(path, params=TMParams(), frame=frame)
            dg = -float(pred.membrane_score)
            n_tm = len(pred.segments or [])
        except Exception as exc:
            n_tm = -1
            print(f"[bench] soluble {acc}: {exc}")
        rows.append({"accession": acc, "pdb_id": pdb_id, "chain": frame.chain_id,
                     "n_residues": len(frame), "dG_transfer_kcal_per_mol": dg,
                     "n_tm_segments": n_tm,
                     "called_membrane": bool(dg is not None and dg <= -10.0),
                     "protein_name": (entry.get("proteinDescription") or {})
                     .get("recommendedName", {}).get("fullName", {}).get("value", "")})
        print(f"[bench] soluble {acc} {pdb_id} dG={dg} nTM={n_tm}")
    (out_dir / "soluble_control.json").write_text(json.dumps(rows, indent=2))
    if rows:
        dgs = [r["dG_transfer_kcal_per_mol"] for r in rows
               if r["dG_transfer_kcal_per_mol"] is not None]
        fp = sum(1 for r in rows if r["called_membrane"])
        summary = {"n": len(rows), "n_called_membrane": fp,
                   "dG_min": min(dgs) if dgs else None,
                   "dG_max": max(dgs) if dgs else None,
                   "dG_median": round(float(np.median(dgs)), 2) if dgs else None}
        (out_dir / "soluble_control_summary.json").write_text(json.dumps(summary, indent=2))
        print("[bench] soluble summary:", summary)
    return rows


def summarise(results, out_dir: Path):
    names = [m for m in list(METHODS) + ["toppred", "topcons_style_consensus"]
             if any(m in r["methods"] for r in results)]
    lines = ["method,n_proteins,n_ref,n_pred,n_match,macro_recall,macro_precision,"
             "macro_f1,residue_recall,median_max_boundary_error,frac_within_5_residues,"
             "mean_runtime_s"]
    summary = []
    for m in names:
        rows = [r["methods"][m] for r in results if m in r["methods"] and "error" not in r["methods"][m]]
        if not rows:
            continue
        def mean(key):
            vals = [r[key] for r in rows if r.get(key) is not None]
            return round(float(np.mean(vals)), 4) if vals else None
        rec = {
            "method": m, "n_proteins": len(rows),
            "n_ref": int(sum(r["n_ref"] for r in rows)),
            "n_pred": int(sum(r["n_pred"] for r in rows)),
            "n_match": int(sum(r["n_match"] for r in rows)),
            "macro_recall": mean("recall"), "macro_precision": mean("precision"),
            "macro_f1": mean("f1"), "residue_recall": mean("residue_recall"),
            "median_max_boundary_error": mean("median_max_boundary"),
            "frac_within_5_residues": mean("within_tol_frac"),
            "mean_runtime_s": mean("seconds"),
        }
        summary.append(rec)
        lines.append(",".join("" if rec[k] is None else str(rec[k]) for k in
                              ("method", "n_proteins", "n_ref", "n_pred", "n_match",
                               "macro_recall", "macro_precision", "macro_f1", "residue_recall",
                               "median_max_boundary_error", "frac_within_5_residues",
                               "mean_runtime_s")))
    (out_dir / "summary.csv").write_text("\n".join(lines) + "\n")
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="bench/results")
    ap.add_argument("--limit", type=int, default=60)
    ap.add_argument("--query", default="(reviewed:true) AND (ft_transmem:*) "
                                       "AND (organism_id:9606)")
    ap.add_argument("--soluble", action="store_true",
                    help="also run the soluble-protein control cohort")
    ap.add_argument("--workers", type=int, default=1,
                    help="number of parallel worker processes")
    args = ap.parse_args()
    res = run(Path(args.out), args.limit, args.query, workers=args.workers)
    print("\n=== summary ===")
    for rec in summarise(res, Path(args.out)):
        print(rec)
    if args.soluble:
        print("\n=== soluble control ===")
        run_soluble(Path(args.out))

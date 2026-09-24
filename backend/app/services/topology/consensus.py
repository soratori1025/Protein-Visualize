"""
Consensus flow, residue level: TM block (UniProt ...) x SS block (DSSP / STRIDE).

    for i in 0 .. N-1:
        flag_tm = residue i lies in a Transmembrane segment of the TM block
        flag_ss = residue i is a helix in the SS block (a strand inside a beta TM segment)
        flag_ss and flag_tm        -> TM_in   (helix inside the bilayer band)
        flag_ss and not flag_tm    -> TM_C    (helix on the cytoplasmic side)
                                      TM_E    (helix on the extracellular / lumenal side)
        not flag_ss                -> not put in the map

The side of a TM_C / TM_E residue comes from the TOPOLOGY array (UniProt
"Topological domain", or the side the orchestrator inferred), never from DSSP: DSSP
only knows hydrogen bonds, not which face of the membrane a residue is on.

The map is RAW: one entry per residue, exactly as the rule says. Two chemistry rules
are applied only when TM_in residues are turned into crossings for the diagram
(``tm_in_runs``):
  * a hole of 1-2 non-helical residues inside a TM (kink, or a DSSP turn at a
    proline) does not break the helix  (MIN_UNWOUND = 3)
  * a TM_in fragment shorter than MIN_FRAGMENT_LEN is not drawn as its own part.

All arrays are indexed by ResidueFrame POSITION (0-based), so position i is the same
residue in the TM array, the SS array and the frame. Different lengths mean the
arrays were built on different residue lists: that is an error, never truncated.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Callable, Optional, Sequence

from app.services.topology import labels as L
from app.services.topology.transitions import MIN_FRAGMENT_LEN, MIN_UNWOUND

TM_IN = "TM_in"
TM_C = "TM_C"
TM_E = "TM_E"

# DSSP/STRIDE helix codes. H = alpha, G = 3-10, I = pi. G and I are kept on purpose:
# TM helices often end in a 3-10 turn and carry pi-bulges in the middle, and DSSP >= 2.1
# gives pi priority over alpha, so an H-only rule would cut a normal TM helix in two
# at every pi-bulge (a false TM1a/1b). The raw code stays in the map (ss_raw).
HELIX_CODES = frozenset("HGI")
# Strands of beta-barrel TM segments. B (isolated bridge) is not a strand.
STRAND_CODES = frozenset("E")


class IndexMismatchError(ValueError):
    """The TM and SS arrays were not built on the same residue list."""


@dataclass(frozen=True)
class ConsensusEntry:
    label: str                   # TM_in / TM_C / TM_E
    ss_raw: str                  # DSSP/STRIDE code of the residue (H G I E)
    tm_segment: Optional[int]    # 0-based TM segment of the TM block (TM_in only)


def segment_ids(n: int, tm_labels: Sequence[str],
                tm_segments: Optional[Sequence[tuple[int, int]]] = None) -> list[int]:
    """TM segment index per position, -1 outside. With ``tm_segments`` (the TM block's
    own features) two TMs that touch (TM1 ends at 40, TM2 starts at 41) stay two
    segments; from the labels alone they would look like one."""
    ids = [-1] * n
    if tm_segments:
        for k, (s, e) in enumerate(sorted(tm_segments)):
            for i in range(max(0, s), min(n - 1, e) + 1):
                ids[i] = k
        return ids
    k = -1
    for i in range(n):
        if L.is_tm(tm_labels[i]):
            if i == 0 or not L.is_tm(tm_labels[i - 1]):
                k += 1
            ids[i] = k
    return ids


def build_consensus_map(tm_labels: Sequence[str], ss_codes: Sequence[Optional[str]], *,
                        sides: Optional[Sequence[Optional[str]]] = None,
                        tm_segments: Optional[Sequence[tuple[int, int]]] = None,
                        strand_segments: frozenset[int] = frozenset(),
                        helix_codes: frozenset[str] = HELIX_CODES,
                        ) -> tuple[dict[int, ConsensusEntry], Counter]:
    """The residue loop. Returns ({position: entry}, counts of skipped cases).

    tm_labels   TM block label per position ("Transmembrane", "Cytoplasmic", ...)
    ss_codes    raw DSSP/STRIDE code per position, None = no assignment
    sides       side per position for TM_C / TM_E ("Cytoplasmic" / "Extracellular");
                default: the side written in tm_labels (UniProt topological domain)
    tm_segments (start, end) of each TM segment of the TM block (keeps touching TMs apart)
    strand_segments  TM segments that are beta strands (barrels): flag_ss uses E there
    """
    n = len(tm_labels)
    if len(ss_codes) != n:
        raise IndexMismatchError(f"TM array has {n} residues, SS array {len(ss_codes)}: "
                                 "both must be built on the same residue frame")
    if sides is not None and len(sides) != n:
        raise IndexMismatchError(f"side array has {len(sides)} residues, expected {n}")
    seg = segment_ids(n, tm_labels, tm_segments)

    cmap: dict[int, ConsensusEntry] = {}
    skipped: Counter = Counter()
    for i in range(n):
        flag_tm = seg[i] >= 0
        code = ss_codes[i]
        accepted = STRAND_CODES if flag_tm and seg[i] in strand_segments else helix_codes
        if code is None or code not in accepted:
            continue                                   # not a helix: not in the map
        if flag_tm:
            cmap[i] = ConsensusEntry(TM_IN, code, seg[i])
            continue
        side = _side(sides[i] if sides is not None else tm_labels[i])
        if side == L.CYTO:
            cmap[i] = ConsensusEntry(TM_C, code, None)
        elif side == L.EXTRA:
            cmap[i] = ConsensusEntry(TM_E, code, None)
        else:
            # helix outside the TM band whose side is unknown: Intramembrane (re-entrant
            # half helix, pore helix), Signal peptide, or no topology information
            base = L.normalize_base_label(tm_labels[i])
            skipped[base if base in (L.INTRA, L.SIGNAL) else "no side"] += 1
    return cmap, skipped


def _side(label: Optional[str]) -> Optional[str]:
    if not label:
        return None
    first = label.split()[0]
    if first in L.KNOWN_SIDES:
        return first
    return L.side_from_text(label) if not L.is_tm(label) and L.INTRA not in label else None


def tm_in_runs(cmap: dict[int, ConsensusEntry], n: int,
               breaks: Optional[Sequence[bool]] = None) -> dict[int, list[tuple[int, int]]]:
    """TM_in positions -> maximal runs per TM segment ({segment: [(start, end), ...]}).
    A run never continues across a chain break (missing residues are not helix)."""
    runs: dict[int, list[tuple[int, int]]] = {}
    i = 0
    while i < n:
        e = cmap.get(i)
        if e is None or e.label != TM_IN:
            i += 1
            continue
        assert e.tm_segment is not None
        j = i + 1
        while (j < n and (f := cmap.get(j)) is not None and f.label == TM_IN
               and f.tm_segment == e.tm_segment and not (breaks is not None and breaks[j])):
            j += 1
        runs.setdefault(e.tm_segment, []).append((i, j - 1))
        i = j
    return runs


def join_kinks(runs: list[tuple[int, int]], breaks: Optional[Sequence[bool]] = None,
               min_unwound: int = MIN_UNWOUND,
               is_hairpin: Optional[Callable[[tuple[int, int], tuple[int, int]], bool]] = None,
               ) -> list[tuple[int, int]]:
    """Runs separated by fewer than ``min_unwound`` residues (and no chain break) are
    one helix with a kink - unless ``is_hairpin(last_run, run)`` says the chain turns
    back (MFS / LeuT hairpins often have only a 1-2 residue turn)."""
    out: list[tuple[int, int]] = []
    last: Optional[tuple[int, int]] = None          # last raw run of the current piece
    for a, b in sorted(runs):
        if out:
            gap = a - out[-1][1] - 1
            broken = breaks is not None and any(breaks[k] for k in range(out[-1][1] + 1, a + 1))
            turns = is_hairpin is not None and last is not None and is_hairpin(last, (a, b))
            if gap < min_unwound and not broken and not turns:
                out[-1] = (out[-1][0], b)
                if b - a + 1 >= MIN_FRAGMENT_LEN:
                    last = (a, b)
                continue
        out.append((a, b))
        last = (a, b)
    return out


def drop_short_fragments(runs: list[tuple[int, int]],
                         min_len: int = MIN_FRAGMENT_LEN) -> list[tuple[int, int]]:
    """Fragments shorter than ``min_len`` are not drawn as parts when a longer one
    exists (they stay TM_in in the raw map)."""
    long_ = [(a, b) for a, b in runs if b - a + 1 >= min_len]
    return long_ if long_ else list(runs)

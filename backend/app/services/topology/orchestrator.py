"""
Topology orchestrator: runs one TM block and one SS block and merges them with one
of three flows.

INDEXING CONTRACT
  The structure is parsed ONCE into a ResidueFrame (residues.py). The TM provider
  returns one label per frame POSITION, the SS provider returns one code per frame
  POSITION (mapped by residue key, verified by residue identity, alignment fallback),
  and every flow below works on those position arrays only. Author residue numbers
  come back only when regions are emitted. Nothing is looked up by a bare int
  residue number and nothing is rebuilt with ``range(start, end + 1)``, so the two
  blocks can never disagree about which residue is which.

FLOWS (each maps TM segments + SS codes -> final TM spans)
  tm_then_ss      TM boundaries are kept exactly; the majority SS class inside each
                  segment names it (Alpha Helix / Beta Strand / Loop).
  ss_then_tm      Each TM segment is snapped to the SS elements of its dominant class
                  that overlap it: trimmed where the ends are coil, extended along the
                  element by at most ``max_snap`` residues and never so far that the
                  loop to the neighbouring segment disappears.
  parallel_merge  Intersection: a segment keeps only the stretch from its first to its
                  last residue of the dominant SS class. A short interior break (a
                  proline kink, a pi-bulge read as 'T') stays inside - one crossing is
                  never split by a kink. A segment with no helix/strand residue is dropped.

  In both SS-driven flows a segment that contains several elements (a slab or a
  hydropathy window over a helical hairpin) is split into one crossing per element
  group: elements stay together only across a short (<= ``params.broken_gap_max``),
  unbroken gap AND when the chain keeps its direction (kink / pi-bulge); antiparallel
  elements are separate crossings however short the turn. Pieces shorter than
  ``min_tm_core`` are not promoted to crossings.

Freed residues and loops without a side are resolved afterwards: the provider's own
side labels win; otherwise sides alternate across each crossing and the orientation
comes from the positive-inside rule.
"""
from __future__ import annotations
from enum import Enum
from pathlib import Path
from typing import Optional, Sequence
import numpy as np
from app.core.constants import MAX_SNAP, MIN_TM_CORE, POSITIVE_RESIDUES
from app.schemas.topology import TMParams, TopologyRegion, TopologyResponse
from app.services.topology import labels as L
from app.services.topology.providers.ss.base import SSProvider
from app.services.topology.providers.tm.base import TMPrediction, TMProvider
from app.services.topology.residues import (
    ResidueFrame, SSAssignment, SSRecord, load_residue_frame, map_ss_records, parse_resnum,
)

# Loop residues within this distance of a TM end vote in the positive-inside rule.
POSITIVE_INSIDE_FLANK = 15
MIN_SS_COVERAGE_WARN = 0.8

Span = tuple[int, int, Optional[str]]      # (start_pos, end_pos, 'H' | 'E' | 'L' | None)


class FlowType(str, Enum):
    TM_THEN_SS = "tm_then_ss"
    SS_THEN_TM = "ss_then_tm"
    PARALLEL_MERGE = "parallel_merge"


class TopologyOrchestrator:
    def __init__(self, tm_provider: TMProvider, ss_provider: Optional[SSProvider],
                 flow_type: str, *, max_snap: int = MAX_SNAP, min_tm_core: int = MIN_TM_CORE,
                 annotate_extramembrane_ss: bool = True):
        try:
            self.flow = FlowType(flow_type)
        except ValueError:
            raise ValueError(f"Unknown flow_type {flow_type!r}; expected one of "
                             f"{[f.value for f in FlowType]}") from None
        self.flow_type = self.flow.value
        self.tm_provider = tm_provider
        self.ss_provider = ss_provider
        self.max_snap = max_snap
        self.min_tm_core = min_tm_core
        self.annotate_extramembrane_ss = annotate_extramembrane_ss

    # ------------------------------------------------------------------ public
    def execute(self, file_path: Path, params: Optional[TMParams] = None,
                chain_id: Optional[str] = None, **kwargs) -> TopologyResponse:
        params = params if params is not None else TMParams()
        file_path = Path(file_path)
        frame = load_residue_frame(file_path, chain_id)
        warns: list[str] = []
        if len(frame) == 0:
            return self._response(frame, params, TMPrediction(), "No amino-acid residues",
                                  self._tm_name(), [], ["no amino-acid residues in the chain"])

        # 1. TM block, on the shared frame
        tm_pred = self.tm_provider.predict_tm(file_path, params=params, frame=frame, **kwargs)
        warns.extend(tm_pred.warnings)
        base, segments, tm_warns = self._tm_on_frame(tm_pred, frame)
        warns.extend(tm_warns)
        if not segments:
            return self._response(frame, params, tm_pred,
                                  f"No transmembrane segments ({tm_pred.labeler or self._tm_name()})",
                                  tm_pred.labeler or self._tm_name(), [], warns)

        # 2. SS block, on the same frame
        coarse, ss_name, ss_warns = self._ss_on_frame(file_path, frame)
        warns.extend(ss_warns)
        labeler = f"{tm_pred.labeler or self._tm_name()} + {ss_name} ({self.flow_type})"

        # 3. flow -> final TM spans (positions)
        # Convert to list to satisfy Sequence[bool] type hints, and since pure-Python indexing is faster
        breaks = frame.chain_breaks().tolist()
        if coarse is None:
            spans: list[Span] = [(s, e, None) for s, e in segments]
        elif self.flow is FlowType.TM_THEN_SS:
            spans = self._flow_tm_then_ss(segments, coarse)
        elif self.flow is FlowType.SS_THEN_TM:
            spans = self._flow_ss_then_tm(segments, coarse, breaks, params.broken_gap_max,
                                          frame.coords)
        else:
            spans, dropped = self._flow_parallel_merge(segments, coarse, breaks,
                                                       params.broken_gap_max, frame.coords)
            if dropped:
                warns.append(f"parallel_merge: {dropped} TM segment(s) had no helix/strand "
                             f"residue in {ss_name} and were dropped")
        if not spans:
            return self._response(frame, params, tm_pred,
                                  "No transmembrane segment supported by secondary structure",
                                  labeler, [], warns)

        # 4. paint spans, resolve sides, annotate loops, emit regions
        labels, groups = self._paint(base, spans)
        labels = self._resolve_sides(labels, spans, frame)
        if coarse is not None and self.annotate_extramembrane_ss:
            labels = L.label_extramembrane_ss(labels, coarse, breaks=breaks)
        regions = L.build_regions(frame, labels, groups)
        return self._response(frame, params, tm_pred, "Predicted Topology", labeler, regions, warns)

    # ------------------------------------------------------- blocks -> frame
    def _tm_name(self) -> str:
        return self.tm_provider.__class__.__name__.replace("Provider", "")

    def _tm_on_frame(self, tm_pred: TMPrediction, frame: ResidueFrame):
        """-> (base labels per position, sorted non-overlapping TM segments, warnings)."""
        n = len(frame)
        warns: list[str] = []
        painted_segments: list[tuple[int, int]] = []
        labels_ok = tm_pred.labels is None or len(tm_pred.labels) == n
        if tm_pred.labels is not None and labels_ok:
            base = [L.normalize_base_label(x) for x in tm_pred.labels]
        else:
            if not labels_ok and tm_pred.labels is not None:   # built on another residue list: positions are meaningless here
                warns.append(f"TM provider returned {len(tm_pred.labels)} labels for {n} "
                             "residues; rebuilt from its regions")
            base, painted_segments, outside = self._paint_regions(tm_pred.regions, frame)
            if outside:
                warns.append(f"{outside} TM-provider region(s) lie outside the residues of "
                             f"chain {frame.chain_id} and were ignored")

        use_given = tm_pred.segments is not None and labels_ok
        raw_segments = tm_pred.segments if tm_pred.segments is not None and labels_ok else painted_segments
        segments = [(max(0, int(s)), min(n - 1, int(e))) for s, e in raw_segments
                    if int(e) >= 0 and int(s) < n and int(s) <= int(e)]
        if not segments:
            segments = self._runs(base, L.TM)
        segments = self._normalize_segments(segments)
        # labels and segments must agree: TM positions outside any segment become one
        covered = [False] * n
        for s, e in segments:
            for k in range(s, e + 1):
                covered[k] = True
                base[k] = L.TM
        extra = [(s, e) for s, e in self._runs(base, L.TM)
                 if not any(covered[k] for k in range(s, e + 1))]
        segments = self._normalize_segments(segments + extra)
        return base, segments, warns

    @staticmethod
    def _paint_regions(regions: Sequence[TopologyRegion], frame: ResidueFrame):
        """Fallback for providers that only give author-numbered regions: paint them
        on OBSERVED positions (never range(start, end + 1)); TM painted last."""
        n = len(frame)
        base = [L.UNASSIGNED] * n
        segments, outside = [], 0

        def label_of(r):
            if r.type == L.TM or (r.type != "Topological domain" and L.TM in (r.description or "")):
                return L.TM
            if r.type == L.INTRA:
                return L.INTRA
            if r.type == L.SIGNAL:
                return L.SIGNAL
            if r.side in L.KNOWN_SIDES:
                return r.side
            return L.side_from_text(r.description) or L.normalize_base_label(r.description)

        ordered = sorted(regions, key=lambda r: (label_of(r) == L.TM, r.start))
        for r in ordered:
            lo, hi = min(r.start, r.end), max(r.start, r.end)
            span = frame.span_to_positions(lo, hi, r.start_icode if r.start <= r.end else None,
                                           r.end_icode if r.start <= r.end else None)
            if span is None:
                outside += 1
                continue
            lab = label_of(r)
            for k in range(span[0], span[1] + 1):
                base[k] = lab
            if lab == L.TM:
                segments.append(span)
        return base, segments, outside

    def _ss_on_frame(self, file_path: Path, frame: ResidueFrame):
        """-> (coarse 'H'/'E'/'C'/None per position, or None if unusable; name; warnings)."""
        if self.ss_provider is None:
            return None, "None", ["no SS provider: TM segments are reported without an SS class"]
        name = getattr(self.ss_provider, "name", None) or \
            self.ss_provider.__class__.__name__.replace("Provider", "")
        try:
            if hasattr(self.ss_provider, "assign"):
                result: SSAssignment = self.ss_provider.assign(file_path, frame)
            else:   # third-party provider with the old Dict[int, str] contract
                legacy = self.ss_provider.get_secondary_structure(file_path)
                records = []
                for key, code in legacy.items():
                    num = parse_resnum(key)
                    if num is not None:
                        records.append(SSRecord(num[0], num[1], None, str(code or "C"), None))
                result = map_ss_records(frame, records, labeler=name)
        except Exception as error:  # noqa: BLE001 - external binary may be missing
            print(f"[topology] {name} failed: {error}")
            return None, f"{name} (unavailable)", [
                f"{name} failed ({error}); TM segments are reported without an SS class"]

        warns = list(result.warnings)
        if result.matched == 0:
            warns.append(f"{name}: no residue could be matched to chain {frame.chain_id}; "
                         "SS ignored")
            return None, f"{name} (unmatched)", warns
        if result.coverage < MIN_SS_COVERAGE_WARN:
            warns.append(f"{name}: only {result.matched}/{len(frame)} residues have an SS "
                         "assignment")
        return [L.coarse_ss(c) for c in result.codes], name, warns

    # ------------------------------------------------------------------ flows
    @staticmethod
    def _majority(codes: Sequence[Optional[str]]) -> str:
        h, e = codes.count("H"), codes.count("E")
        if h == 0 and e == 0:
            return "L"
        return "H" if h >= e else "E"

    def _flow_tm_then_ss(self, segments, coarse) -> list[Span]:
        return [(s, e, self._majority(coarse[s:e + 1])) for s, e in segments]

    def _flow_ss_then_tm(self, segments, coarse, breaks, gap_max, coords=None) -> list[Span]:
        n = len(coarse)
        elements = self._ss_elements(coarse, breaks)
        spans: list[Span] = []
        for idx, (s0, e0) in enumerate(segments):
            overlapping = [(a, b, c) for a, b, c in elements if a <= e0 and b >= s0]
            if not overlapping:
                spans.append((s0, e0, "L"))
                continue
            h = sum(min(b, e0) - max(a, s0) + 1 for a, b, c in overlapping if c == "H")
            st = sum(min(b, e0) - max(a, s0) + 1 for a, b, c in overlapping if c == "E")
            cls = "H" if h >= st else "E"
            mine = [(a, b) for a, b, c in overlapping if c == cls]

            prev_end = segments[idx - 1][1] if idx > 0 else -1
            next_start = segments[idx + 1][0] if idx + 1 < len(segments) else n
            left_gap, right_gap = s0 - prev_end - 1, next_start - e0 - 1
            left_budget = max(0, min(self.max_snap,
                                     left_gap // 2 if idx == 0 else (left_gap - 1) // 2))
            right_budget = max(0, min(self.max_snap, right_gap // 2
                                      if idx + 1 == len(segments) else (right_gap - 1) // 2))
            lo_cap, hi_cap = s0 - left_budget, e0 + right_budget

            # One TM segment may hold several elements (a coarse hydropathy call over two
            # helices and their short loop): elements separated by more than `gap_max`
            # coil residues become separate crossings; a shorter break (kink) does not.
            pieces = [(max(a, lo_cap), min(b, hi_cap))
                      for a, b in self._group_positions(mine, gap_max, breaks, coords)]
            keep_len = min(self.min_tm_core, e0 - s0 + 1)
            pieces = [(a, b) for a, b in pieces if b - a + 1 >= keep_len]
            if not pieces:
                spans.append((s0, e0, cls))       # SS support too thin to move the bounds
            else:
                spans.extend((a, b, cls) for a, b in pieces)
        return spans

    def _flow_parallel_merge(self, segments, coarse, breaks, gap_max,
                             coords=None) -> tuple[list[Span], int]:
        spans: list[Span] = []
        dropped = 0
        for s0, e0 in segments:
            cls = self._majority(coarse[s0:e0 + 1])
            if cls == "L":
                dropped += 1
                continue
            support = [(a, b) for a, b, c in self._ss_elements(coarse[:e0 + 1], breaks)
                       if c == cls and b >= s0]
            support = [(max(a, s0), b) for a, b in support]
            pieces = self._group_positions(support, gap_max, breaks, coords)
            keep_len = min(self.min_tm_core, e0 - s0 + 1)
            kept = [p for p in pieces if p[1] - p[0] + 1 >= keep_len]
            if not kept:
                kept = [max(pieces, key=lambda p: p[1] - p[0])]
            spans.extend((a, b, cls) for a, b in kept)
        return spans, dropped

    @staticmethod
    def _group_positions(intervals, gap_max, breaks, coords=None) -> list[tuple[int, int]]:
        """Group sorted inclusive intervals (SS elements, or single residues) into
        crossings. The next interval joins the current group only if the gap between
        them is <= gap_max residues, contains no chain break, and - when both pieces are
        long enough to have a direction - the chain keeps going the same way (a kink or
        a pi-bulge). Two antiparallel elements are a hairpin, i.e. two crossings, even
        when the turn between them is only 1-2 residues (common in MFS transporters)."""
        groups: list[tuple[int, int]] = []
        last: Optional[tuple[int, int]] = None           # last element of current group

        def direction(a, b):
            if coords is None or b - a < 3:
                return None
            v = coords[b] - coords[a]
            norm = float(np.linalg.norm(v))
            return v / norm if np.isfinite(norm) and norm > 0 else None

        for a, b in sorted(intervals):
            if groups:
                ga, gb = groups[-1]
                gap = a - gb - 1
                joinable = gap <= gap_max and not any(breaks[k] for k in range(gb + 1, a + 1))
                if joinable and last is not None:
                    u, v = direction(*last), direction(a, b)
                    if u is not None and v is not None and float(u @ v) < 0.0:
                        joinable = False                    # antiparallel -> hairpin
                if joinable:
                    groups[-1] = (ga, max(gb, b))
                    if b - a >= 3:                          # keep a piece with a direction
                        last = (a, b)
                    continue
            groups.append((a, b))
            last = (a, b)
        return groups

    # ---------------------------------------------------------------- helpers
    @staticmethod
    def _runs(labels: Sequence[str], value: str) -> list[tuple[int, int]]:
        runs, i, n = [], 0, len(labels)
        while i < n:
            if labels[i] != value:
                i += 1
                continue
            j = i
            while j < n and labels[j] == value:
                j += 1
            runs.append((i, j - 1))
            i = j
        return runs

    @staticmethod
    def _normalize_segments(segments) -> list[tuple[int, int]]:
        out: list[tuple[int, int]] = []
        for s, e in sorted(set(segments)):
            if out and s <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], e))    # overlap -> one segment
            else:
                out.append((s, e))
        return out

    @staticmethod
    def _ss_elements(coarse, breaks) -> list[tuple[int, int, str]]:
        """Maximal helix / strand runs; an element never continues across a chain break."""
        elements, i, n = [], 0, len(coarse)
        while i < n:
            c = coarse[i]
            if c not in ("H", "E"):
                i += 1
                continue
            j = i + 1
            while j < n and coarse[j] == c and not breaks[j]:
                j += 1
            elements.append((i, j - 1, c))
            i = j
        return elements

    @staticmethod
    def _paint(base: Sequence[str], spans: Sequence[Span]):
        """TM positions not covered by a final span are freed (Unassigned, side decided
        later); span positions get the span label and a unique group id."""
        labels = [L.UNASSIGNED if lab == L.TM else lab for lab in base]
        groups = [-1] * len(base)
        for g, (s, e, cls) in enumerate(spans):
            for k in range(s, e + 1):
                labels[k] = L.TM_CLASS_LABEL[cls]
                groups[k] = g
        return labels, groups

    @staticmethod
    def _opposite(side: str) -> str:
        return L.EXTRA if side == L.CYTO else L.CYTO

    def _resolve_sides(self, labels: list[str], spans: Sequence[Span],
                       frame: ResidueFrame) -> list[str]:
        labels = L.apply_positive_inside_rule(labels, frame.names)
        n = len(labels)
        loops = []
        i = 0
        while i < n:
            if L.is_tm(labels[i]):
                i += 1
                continue
            j = i
            while j < n and not L.is_tm(labels[j]):
                j += 1
            loops.append((i, j - 1))
            i = j
        if not loops:
            return labels
        crossings_before = [sum(1 for s, _, _ in spans if s < a) for a, _ in loops]

        entry: list[Optional[str]] = []
        exit_: list[Optional[str]] = []
        for a, b in loops:
            known = [p for p in range(a, b + 1) if labels[p] in L.KNOWN_SIDES]
            if not known:
                entry.append(None)
                exit_.append(None)
                continue
            entry.append(labels[known[0]])
            exit_.append(labels[known[-1]])
            for p in range(a, b + 1):
                if labels[p] == L.UNASSIGNED:
                    nearest = min(known, key=lambda q: (abs(q - p), q > p))
                    labels[p] = labels[nearest]

        anchors = [k for k, side in enumerate(entry) if side is not None]
        if anchors:
            inferred = {}
            for k in range(len(loops)):
                if entry[k] is not None:
                    continue
                a = min(anchors, key=lambda x: (abs(x - k), x > k))
                ref = exit_[a] if a < k else entry[a]
                if ref is not None:
                    flips = abs(crossings_before[k] - crossings_before[a])
                    inferred[k] = ref if flips % 2 == 0 else self._opposite(ref)
        else:
            inferred = self._positive_inside_by_parity(loops, crossings_before, frame, labels)

        for k, side in inferred.items():
            a, b = loops[k]
            for p in range(a, b + 1):
                if labels[p] == L.UNASSIGNED:
                    labels[p] = side
        return labels

    @staticmethod
    def _positive_inside_by_parity(loops, crossings_before, frame, labels) -> dict:
        """No side known anywhere (sequence-only TM calls): loops alternate across
        crossings; the parity whose membrane-flanking residues are richer in Lys/Arg is
        cytoplasmic (von Heijne). Ties -> N-terminal loop parity inside."""
        n = len(labels)
        pos = {0: 0, 1: 0}
        tot = {0: 0, 1: 0}
        names = frame.names
        for k, (a, b) in enumerate(loops):
            parity = crossings_before[k] % 2
            idx = set()
            if a > 0:                                   # TM span before the loop
                idx.update(range(a, min(b, a + POSITIVE_INSIDE_FLANK - 1) + 1))
            if b < n - 1:                               # TM span after the loop
                idx.update(range(max(a, b - POSITIVE_INSIDE_FLANK + 1), b + 1))
            pos[parity] += sum(1 for p in idx if names[p] in POSITIVE_RESIDUES)
            tot[parity] += len(idx)
        dens = {p: (pos[p] / tot[p] if tot[p] else 0.0) for p in (0, 1)}
        inside = 0 if dens[0] >= dens[1] else 1
        return {k: (L.CYTO if crossings_before[k] % 2 == inside else L.EXTRA)
                for k in range(len(loops))}

    def _response(self, frame, params, tm_pred, name, labeler, regions, warns) -> TopologyResponse:
        normal = tm_pred.membrane_normal
        return TopologyResponse(
            uniprot_id="CALCULATED",
            protein_name=name,
            gene_name="",
            organism="Computed",
            membrane_score=round(float(tm_pred.membrane_score or 0.0), 3),
            membrane_normal=[round(float(x), 4) for x in normal] if normal else None,
            labeler=labeler,
            parameters_used={**params.to_response_dict(), "flow": self.flow_type,
                             "max_snap": self.max_snap},
            regions=regions,
            chain_id=frame.chain_id,
            warnings=list(dict.fromkeys(warns)),
        )

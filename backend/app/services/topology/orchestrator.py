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

LAYERS (all per residue, on the frame)
  TM evidence        the TM block (slab geometry, Kyte-Doolittle, UniProt ...)
  SS evidence        DSSP / STRIDE, raw 8-state codes kept, coarse H/E/C for rules
  Membrane geometry  where the bilayer is (membrane.py): file planes (OPM/PPM DUM
                     atoms) > axes of the TM segments; depth + CORE/EDGE/OUT zone.
  A DSSP boundary is never used as a membrane boundary: SS only describes the part
  of a TM candidate that the geometry places in the bilayer.

FLOWS (API value / alias)
  tm_then_ss / tm_first
      TM boundaries are kept. Each segment is classified from its SS content:
      Alpha Helix / Beta Strand when one class dominates, Irregular when mixed or
      unassigned (e.g. H + E elements, or coil only).
  ss_then_tm / structure_guided
      Bounded refinement. Candidate window = segment +- min(max_snap, half the loop
      to the neighbour). SS elements inside the window are grouped into crossings:
      elements stay together only across a short (<= broken_gap_max), unbroken gap
      that stays inside the membrane AND when the chain keeps its direction; each
      group is clipped to the membrane envelope. A group that crosses (span >=
      min_cross_span_frac x thickness) is a TM crossing; one that dips into the core
      without crossing is Intramembrane (re-entrant / half helix). A segment with no
      SS support is kept as Irregular if the geometry puts it in the membrane,
      otherwise dropped with a warning. Without coordinates: SS-only snapping.
  parallel_merge / consensus
      Residue-level agreement: TM evidence AND SS (dominant class) AND membrane
      envelope. Kinks inside a crossing are tolerated, hairpins split. Disagreement
      lowers confidence instead of deleting: a segment with no SS support stays as
      Irregular (low confidence) unless the membrane geometry rejects it.

After the flow: sides (provider first, else alternation + positive-inside),
extramembrane SS, interfacial helices (helices lying in the interface band parallel
to the membrane), per-residue evidence matrix + confidence, domain type and a
topology validation that only warns.

Freed residues and loops without a side are resolved afterwards: the provider's own
side labels win; otherwise sides alternate across each crossing and the orientation
comes from the positive-inside rule.
"""
from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np

from app.core.constants import MAX_SNAP, MIN_TM_CORE, POSITIVE_RESIDUES
from app.schemas.topology import ResidueAnnotation, TMParams, TopologyRegion, TopologyResponse
from app.services.topology import labels as L
from app.services.topology.membrane import EDGE_WIDTH, MembraneFrame, build_membrane
from app.services.topology.transitions import (
    MERGING, MIN_FRAGMENT_LEN, MIN_UNWOUND, TransitionEvidence, classify_membrane_transition,
)
from app.services.topology.providers.ss.base import SSProvider
from app.services.topology.providers.tm.base import TMPrediction, TMProvider
from app.services.topology.residues import (
    ResidueFrame, SSAssignment, SSRecord, load_residue_frame, map_ss_records, parse_resnum,
)

# Loop residues within this distance of a TM end vote in the positive-inside rule.
POSITIVE_INSIDE_FLANK = 15
MIN_SS_COVERAGE_WARN = 0.8

Span = tuple[int, int, Optional[str]]      # (start_pos, end_pos, 'H' | 'E' | 'I' | None)


class FlowType(str, Enum):
    TM_THEN_SS = "tm_then_ss"          # a.k.a. tm_first
    SS_THEN_TM = "ss_then_tm"          # a.k.a. structure_guided
    PARALLEL_MERGE = "parallel_merge"  # a.k.a. consensus


FLOW_ALIASES = {"tm_first": "tm_then_ss", "structure_guided": "ss_then_tm",
                "consensus": "parallel_merge"}

# Plausible lengths of one crossing (residues) - validation warnings only.
TM_LENGTH_RANGE = {"H": (14, 40), "E": (5, 16)}
# Broken crossings (TM1a/1b ...) are decided by transitions.classify_membrane_transition.
INTERFACIAL_MIN_LEN = 6


class TopologyOrchestrator:
    def __init__(self, tm_provider: TMProvider, ss_provider: Optional[SSProvider],
                 flow_type: str, *, max_snap: int = MAX_SNAP, min_tm_core: int = MIN_TM_CORE,
                 annotate_extramembrane_ss: bool = True, include_residues: bool = True):
        try:
            self.flow = FlowType(FLOW_ALIASES.get(flow_type, flow_type))
        except ValueError:
            raise ValueError(f"Unknown flow_type {flow_type!r}; expected one of "
                             f"{[f.value for f in FlowType] + list(FLOW_ALIASES)}") from None
        self.flow_type = self.flow.value
        self.tm_provider = tm_provider
        self.ss_provider = ss_provider
        self.max_snap = max_snap
        self.min_tm_core = min_tm_core
        self.annotate_extramembrane_ss = annotate_extramembrane_ss
        self.include_residues = include_residues

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

        # 2. SS block and membrane geometry, on the same frame
        coarse, raw_ss, ss_name, ss_warns = self._ss_on_frame(file_path, frame)
        warns.extend(ss_warns)
        membrane = build_membrane(frame, segments, params.membrane_thickness / 2.0)
        if membrane is None:
            warns.append("no membrane geometry (coordinates missing): boundaries from TM/SS only")
        elif coarse is not None:
            # centre the estimated bilayer on the helices/strands that the TM block
            # proposed (long SS elements overlapping a TM segment)
            spanning = [(a, b) for a, b, _ in self._ss_elements(coarse, frame.chain_breaks())
                        if b - a + 1 >= 8 and any(a <= e and b >= s for s, e in segments)
                        and membrane.span(a, b) >= params.min_cross_span_frac * 2 * membrane.half_thickness]
            membrane = membrane.recentered(spanning)
        labeler = f"{tm_pred.labeler or self._tm_name()} + {ss_name} ({self.flow_type})"

        # 3. flow -> final TM spans (+ intramembrane pieces)
        breaks = frame.chain_breaks()
        intramembrane: list[Span] = []
        if coarse is None:
            spans: list[Span] = [(s, e, None) for s, e in segments]
        elif self.flow is FlowType.TM_THEN_SS:
            spans = self._flow_tm_then_ss(segments, coarse)
            if membrane is not None:
                spans, notes = self._flag_multi_crossing(spans, coarse, breaks, frame,
                                                         membrane, params)
                warns.extend(notes)
        elif self.flow is FlowType.SS_THEN_TM:
            if membrane is not None:
                spans, intramembrane, notes = self._flow_structure_guided(
                    segments, coarse, breaks, frame.coords, membrane, params)
                warns.extend(notes)
            else:
                spans = self._flow_ss_then_tm(segments, coarse, breaks, params.broken_gap_max,
                                              frame.coords)
        else:
            spans, intramembrane, kept_irregular, dropped = self._flow_parallel_merge(
                segments, coarse, breaks, params.broken_gap_max, frame.coords, membrane, params)
            if kept_irregular:
                warns.append(f"consensus: {kept_irregular} TM segment(s) have no helix/strand "
                             f"support in {ss_name}; kept as Transmembrane Irregular "
                             "(confidence reduced, not deleted)")
            if dropped:
                warns.append(f"consensus: {dropped} TM segment(s) have neither SS support nor "
                             "a place in the membrane geometry and were dropped")
        if not spans:
            return self._response(frame, params, tm_pred,
                                  "No transmembrane segment supported by the evidence",
                                  labeler, [], warns, membrane=membrane)

        # 4. paint, sides, loop annotation
        spans = sorted(spans)
        broken, transitions = self._find_discontinuities(spans, coarse, breaks, frame,
                                                         membrane, params)
        for idx, evs in transitions.items():
            if any(ev["classification"] == "AMBIGUOUS" for ev in evs):
                warns.append(f"TM{idx + 1}: one fragment crosses the bilayer alone and the next "
                             "only partly - drawn as one continuous crossing, confidence lowered")
        labels, groups = self._paint(base, spans)
        for _, (u0, u1), _ in broken.values():
            for k in range(u0, u1 + 1):
                labels[k] = L.UNWOUND
        for s0, e0, _ in intramembrane:
            for k in range(s0, e0 + 1):
                if not L.is_tm(labels[k]):
                    labels[k] = L.INTRA
        labels = self._resolve_sides(labels, spans, frame)
        if coarse is not None and membrane is not None:
            labels = self._mark_interfacial(labels, coarse, breaks, frame.coords, membrane)
        if coarse is not None and self.annotate_extramembrane_ss:
            labels = L.label_extramembrane_ss(labels, coarse, breaks=list(breaks))

        # 5. evidence matrix, confidence, regions, validation
        conf = self._confidence(labels, spans, base, coarse, membrane, frame)
        regions = L.build_regions(frame, labels, groups)
        self._attach_region_confidence(regions, frame, conf)
        self._attach_crossings(regions, frame, groups, broken, transitions)
        self._attach_topology_labels(regions)
        warns.extend(self._validate(labels, spans, frame))
        residues = (self._residue_matrix(frame, base, raw_ss, coarse, membrane, labels, conf)
                    if self.include_residues else None)
        return self._response(frame, params, tm_pred, "Predicted Topology", labeler, regions,
                              warns, membrane=membrane, domain_type=self._domain_type(spans),
                              residues=residues)

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
            if not labels_ok:   # built on another residue list: positions are meaningless here
                warns.append(f"TM provider returned {len(tm_pred.labels or [])} labels for {n} "
                             "residues; rebuilt from its regions")
            base, painted_segments, outside = self._paint_regions(tm_pred.regions, frame)
            if outside:
                warns.append(f"{outside} TM-provider region(s) lie outside the residues of "
                             f"chain {frame.chain_id} and were ignored")

        use_given = tm_pred.segments is not None and labels_ok
        raw_segments = tm_pred.segments if tm_pred.segments is not None and use_given else painted_segments
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
        """-> (coarse 'H'/'E'/'C'/None per position or None, raw 8-state codes or None,
        name, warnings)."""
        if self.ss_provider is None:
            return None, None, "None", ["no SS provider: TM segments are reported without an SS class"]
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
            return None, None, f"{name} (unavailable)", [
                f"{name} failed ({error}); TM segments are reported without an SS class"]

        warns = list(result.warnings)
        if result.matched == 0:
            warns.append(f"{name}: no residue could be matched to chain {frame.chain_id}; "
                         "SS ignored")
            return None, None, f"{name} (unmatched)", warns
        if result.coverage < MIN_SS_COVERAGE_WARN:
            warns.append(f"{name}: only {result.matched}/{len(frame)} residues have an SS "
                         "assignment")
        return [L.coarse_ss(c) for c in result.codes], list(result.codes), name, warns

    # ------------------------------------------------------------------ flows
    @staticmethod
    def _majority(codes: Sequence[Optional[str]]) -> str:
        """Dominant regular class; 'I' when there is no helix/strand at all."""
        h, e = codes.count("H"), codes.count("E")
        if h == 0 and e == 0:
            return "I"
        return "H" if h >= e else "E"

    @staticmethod
    def _classify_segment(codes: Sequence[Optional[str]]) -> str:
        """tm_first label for a whole segment: a class must cover at least half of the
        segment and the other class less than a quarter; otherwise Irregular
        (mixed H + E elements, or mostly coil)."""
        n = max(1, len(codes))
        h, e = codes.count("H") / n, codes.count("E") / n
        if h >= 0.5 and e < 0.25:
            return "H"
        if e >= 0.5 and h < 0.25:
            return "E"
        return "I"

    def _flow_tm_then_ss(self, segments, coarse) -> list[Span]:
        return [(s, e, self._classify_segment(coarse[s:e + 1])) for s, e in segments]

    def _flow_ss_then_tm(self, segments, coarse, breaks, gap_max, coords=None) -> list[Span]:
        n = len(coarse)
        elements = self._ss_elements(coarse, breaks)
        spans: list[Span] = []
        for idx, (s0, e0) in enumerate(segments):
            overlapping = [(a, b, c) for a, b, c in elements if a <= e0 and b >= s0]
            if not overlapping:
                spans.append((s0, e0, "I"))
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

    def _flow_parallel_merge(self, segments, coarse, breaks, gap_max, coords=None,
                             membrane: Optional[MembraneFrame] = None,
                             params: Optional[TMParams] = None):
        """Consensus: residues where the TM block, the dominant SS class and (when
        available) the membrane envelope agree.
        Returns (spans, intramembrane, kept_irregular, dropped)."""
        spans: list[Span] = []
        intra: list[Span] = []
        kept_irregular = dropped = 0
        params = params if params is not None else TMParams()
        for s0, e0 in segments:
            in_env = (lambda k: True) if membrane is None else membrane.in_envelope
            cls = self._majority(coarse[s0:e0 + 1])
            support = []
            if cls != "I":
                support = [(a, b) for a, b, c in self._ss_elements(coarse[:e0 + 1], breaks)
                           if c == cls and b >= s0]
                support = [(max(a, s0), b) for a, b in support]
                # keep only the part of each element inside the membrane envelope
                clipped = []
                for a, b in support:
                    inside = [k for k in range(a, b + 1) if in_env(k)]
                    if inside:
                        clipped.append((inside[0], inside[-1]))
                support = clipped
            if not support:
                # disagreement: SS says no helix/strand here. Keep the TM call as
                # Irregular unless the membrane geometry also rejects it.
                inside = [k for k in range(s0, e0 + 1) if in_env(k)]
                if len(inside) >= max(1, (e0 - s0 + 1) // 2):
                    spans.append((inside[0], inside[-1], "I"))
                    kept_irregular += 1
                else:
                    dropped += 1
                continue
            if membrane is not None:
                # same crossing logic as structure_guided (direction ACROSS the membrane,
                # so the two halves of a broken helix are never split into two TMs)
                thickness = 2.0 * membrane.half_thickness
                groups = self._group_elements([(a, b, cls) for a, b in support], coarse, breaks,
                                              coords, membrane, gap_max,
                                              params.full_cross_frac * thickness)
                found = False
                for a, b in groups:
                    if membrane.span(a, b) >= params.min_cross_span_frac * thickness:
                        spans.append((a, b, cls))
                        found = True
                    elif np.nanmin(np.abs(membrane.depth[a:b + 1])) <= membrane.half_thickness - EDGE_WIDTH:
                        intra.append((a, b, cls))
                if not found:
                    spans.append((support[0][0], support[-1][1], "I"))
                    kept_irregular += 1
                continue
            pieces = self._group_positions(support, gap_max, breaks, coords)
            keep_len = min(self.min_tm_core, e0 - s0 + 1)
            kept = [p for p in pieces if p[1] - p[0] + 1 >= keep_len]
            if not kept:
                kept = [max(pieces, key=lambda p: p[1] - p[0])]
            spans.extend((a, b, cls) for a, b in kept)
        intra = [p for p in intra if not any(p[0] <= c[1] and p[1] >= c[0] for c in spans)]
        return spans, intra, kept_irregular, dropped

    # ------------------------------------------------ discontinuous crossings
    def _find_discontinuities(self, spans, coarse, breaks, frame, membrane, params):
        """Per crossing: the transitions between its SS fragments (evidence), and for
        BROKEN_TM crossings the parts ((a-part), (unwound), (b-part)). The crossing
        stays ONE crossing; parts only tell the diagram to draw it broken.
        Returns (broken, transitions): {span index: parts}, {span index: [evidence]}."""
        broken, transitions = {}, {}
        if membrane is None or coarse is None:
            return broken, transitions
        thickness = 2.0 * membrane.half_thickness
        elements = self._ss_elements(coarse, breaks)
        for idx, (s0, e0, cls) in enumerate(spans):
            if cls not in ("H", "E"):
                continue
            els = [(max(a, s0), min(b, e0), cls) for a, b, c in elements
                   if c == cls and a <= e0 and b >= s0]
            if len(els) < 2:
                continue
            groups, evidence = self._element_groups(
                els, breaks, frame.coords, membrane,
                max(params.broken_gap_max, MIN_UNWOUND), params.full_cross_frac * thickness,
                params.min_cross_span_frac * thickness, with_evidence=True)
            if len(groups) != 1:
                continue                      # separate crossings (tm_first flags those)
            members, evs = groups[0], evidence[0]
            if evs:
                transitions[idx] = [ev.as_dict() for ev in evs]
            cut = [(ev.gap_residues, i) for i, ev in enumerate(evs) if ev.classification == "BROKEN_TM"]
            if not cut:
                continue
            _, i = max(cut)
            left_end, right_start = members[i][1], members[i + 1][0]
            broken[idx] = ((s0, left_end), (left_end + 1, right_start - 1), (right_start, e0))
        return broken, transitions

    def _flag_multi_crossing(self, spans, coarse, breaks, frame, membrane, params):
        """tm_first keeps the boundaries, but a segment that holds several separate
        crossings (a hairpin under one hydropathy window / slab run) is not "one
        alpha helix": relabel it Irregular and say how many crossings it contains."""
        thickness = 2.0 * membrane.half_thickness
        min_cross = params.min_cross_span_frac * thickness
        out, notes = [], []
        elements = self._ss_elements(coarse, breaks)
        for i, (s0, e0, cls) in enumerate(spans, start=1):
            cand = [(max(a, s0), min(b, e0), c) for a, b, c in elements if a <= e0 and b >= s0]
            groups = self._group_elements(cand, coarse, breaks, frame.coords, membrane,
                                          params.broken_gap_max, params.full_cross_frac * thickness)
            n_cross = sum(1 for a, b in groups if membrane.span(a, b) >= min_cross)
            if n_cross >= 2:
                cls = "I"
                notes.append(f"tm_first: TM{i} ({frame.residues[s0].label}-{frame.residues[e0].label}) "
                             f"contains {n_cross} "
                             "separate crossings - labelled Irregular; structure_guided splits it")
            out.append((s0, e0, cls))
        return out, notes

    # ------------------------------------------------ structure-guided refinement
    def _flow_structure_guided(self, segments, coarse, breaks, coords,
                               membrane: MembraneFrame, params: TMParams):
        """Bounded refinement: SS elements near each TM candidate are grouped into
        crossings; boundaries come from the membrane envelope, never from where a
        helix happens to end. Returns (crossings, intramembrane pieces, notes)."""
        n = len(coarse)
        thickness = 2.0 * membrane.half_thickness
        min_cross = params.min_cross_span_frac * thickness
        full_cross = params.full_cross_frac * thickness
        elements = self._ss_elements(coarse, breaks)
        crossings: list[Span] = []
        intra: list[Span] = []
        notes: list[str] = []
        unsupported = rejected = 0

        for idx, (s0, e0) in enumerate(segments):
            prev_end = segments[idx - 1][1] if idx > 0 else -1
            next_start = segments[idx + 1][0] if idx + 1 < len(segments) else n
            left_gap, right_gap = s0 - prev_end - 1, next_start - e0 - 1
            lo = s0 - max(0, min(self.max_snap, left_gap // 2 if idx == 0 else (left_gap - 1) // 2))
            hi = e0 + max(0, min(self.max_snap, right_gap // 2 if idx + 1 == len(segments)
                                 else (right_gap - 1) // 2))
            cand = [(max(a, lo), min(b, hi), c) for a, b, c in elements if a <= hi and b >= lo]
            groups = self._group_elements(cand, coarse, breaks, coords, membrane,
                                          params.broken_gap_max, full_cross, min_cross)
            found = False
            for a, b in groups:
                inside = [k for k in range(a, b + 1) if membrane.in_envelope(k)]
                if len(inside) < params.min_tm_element_in_slab:
                    continue                    # element does not reach the bilayer
                a2, b2 = inside[0], inside[-1]
                cls = self._majority(coarse[a2:b2 + 1])
                if membrane.span(a2, b2) >= min_cross:
                    crossings.append((a2, b2, cls))
                    found = True
                elif np.nanmin(np.abs(membrane.depth[a2:b2 + 1])) <= membrane.half_thickness - EDGE_WIDTH:
                    intra.append((a2, b2, cls))  # dips into the core without crossing
            if not found:
                inside = [k for k in range(s0, e0 + 1) if membrane.in_envelope(k)]
                if inside and membrane.span(inside[0], inside[-1]) >= min_cross:
                    crossings.append((inside[0], inside[-1], "I"))
                    unsupported += 1
                else:
                    rejected += 1
        # a window can only produce disjoint pieces; drop intramembrane pieces that
        # overlap a crossing from a neighbouring window
        intra = [p for p in intra if not any(p[0] <= c[1] and p[1] >= c[0] for c in crossings)]
        if unsupported:
            notes.append(f"structure_guided: {unsupported} TM segment(s) cross the membrane but "
                         "have no regular helix/strand - labelled Transmembrane Irregular")
        if rejected:
            notes.append(f"structure_guided: {rejected} TM segment(s) do not cross the membrane "
                         f"({membrane.source}) and were not kept as crossings")
        return sorted(crossings), sorted(intra), notes

    def _group_elements(self, cand, coarse, breaks, coords, membrane: MembraneFrame,
                        gap_max: int, full_cross: float,
                        min_cross: Optional[float] = None) -> list[tuple[int, int]]:
        return [(g[0][0], g[-1][1]) for g in
                self._element_groups(cand, breaks, coords, membrane, gap_max, full_cross, min_cross)]

    def _element_groups(
        self,
        cand,
        breaks,
        coords,
        membrane: MembraneFrame,
        gap_max: int,
        full_cross: float,
        min_cross: Optional[float] = None,
        with_evidence: bool = False,
    ) -> Any:
        if min_cross is None:
            min_cross = 0.45 * 2.0 * membrane.half_thickness
        els = sorted((a, b, c) for a, b, c in cand)
        if any(b - a + 1 >= MIN_FRAGMENT_LEN for a, b, _ in els):
            els = [
                (a, b, c)
                for a, b, c in els
                if b - a + 1 >= MIN_FRAGMENT_LEN
            ]
        groups: list[list[tuple[int, int, str]]] = []
        evidence: list[list[TransitionEvidence]] = []
        for a, b, cls in els:
            if groups:
                g = groups[-1]
                prev_a, prev_b, prev_cls = g[-1]
                gap = a - prev_b - 1
                if gap <= gap_max:
                    ev = classify_membrane_transition(
                        (prev_a, prev_b),
                        (a, b),
                        membrane,
                        coords,
                        breaks,
                        min_cross,
                        full_cross,
                        ss_class_a=prev_cls,
                        ss_class_b=cls,
                    )
                    if ev.classification in MERGING:
                        g.append((a, b, cls))
                        evidence[-1].append(ev)
                        continue
            groups.append([(a, b, cls)])
            evidence.append([])
        return (groups, evidence) if with_evidence else groups

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
                flips = abs(crossings_before[k] - crossings_before[a])
                if ref is not None:
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

    # --------------------------------------------------- post-processing layers
    def _mark_interfacial(self, labels, coarse, breaks, coords, membrane: MembraneFrame):
        """Extramembrane helices that lie in the interface band, roughly parallel to the
        membrane plane (amphipathic / interfacial helices, e.g. GPCR helix 8)."""
        out = list(labels)
        n = len(labels)
        i = 0
        while i < n:
            if not (labels[i] in L.KNOWN_SIDES and coarse[i] == "H"):
                i += 1
                continue
            j = i + 1
            while j < n and labels[j] == labels[i] and coarse[j] == "H" and not breaks[j]:
                j += 1
            if j - i >= INTERFACIAL_MIN_LEN:
                zones = [membrane.zone(k) for k in range(i, j)]
                v = coords[j - 1] - coords[i]
                norm = float(np.linalg.norm(v))
                parallel = norm > 0 and abs(float(v @ membrane.normal)) / norm < 0.5
                if parallel and zones.count("EDGE") >= 0.6 * (j - i):
                    for k in range(i, j):
                        out[k] = f"{labels[i]} {L.INTERFACIAL}"
            i = j
        return out

    @staticmethod
    def _level(score: float) -> str:
        return "high" if score >= 0.85 else "medium" if score >= 0.6 else "low"

    def _confidence(self, labels, spans, base, coarse, membrane, frame) -> list[Optional[float]]:
        """Per-residue agreement score in [0, 1] across the available layers."""
        n = len(labels)
        span_class = [None] * n
        for s0, e0, cls in spans:
            for k in range(s0, e0 + 1):
                span_class[k] = cls
        scores: list[Optional[float]] = []
        for k in range(n):
            votes = []
            is_tm = L.is_tm(labels[k])
            votes.append(1.0 if (base[k] == L.TM) == is_tm else 0.0)
            if coarse is not None and coarse[k] is not None:
                if is_tm:
                    want = span_class[k]
                    votes.append(1.0 if want in ("H", "E") and coarse[k] == want else
                                 0.5 if want in ("H", "E") and coarse[k] in ("H", "E") else 0.0)
            if membrane is not None and np.isfinite(membrane.depth[k]):
                inside = membrane.in_envelope(k)
                votes.append(1.0 if inside == is_tm or (not is_tm and membrane.zone(k) != "CORE") else 0.0)
            score = sum(votes) / len(votes)
            if frame.predicted_model and frame.b_factors is not None and frame.b_factors[k] < 50:
                score *= 0.7                     # very low pLDDT: the geometry itself is uncertain
            scores.append(score)
        return scores

    def _attach_region_confidence(self, regions, frame, scores) -> None:
        pos = 0
        for region in regions:
            a = frame.position(region.start, region.start_icode or "")
            b = frame.position(region.end, region.end_icode or "")
            if a is None or b is None:
                continue
            a = max(a, pos)
            vals = [scores[k] for k in range(a, b + 1) if scores[k] is not None]
            if vals:
                region.confidence = self._level(float(np.mean(vals)))
            pos = b + 1

    def _attach_crossings(self, regions, frame, groups, broken, transitions) -> None:
        """crossing = 1-based crossing number; part = 'a'/'b' for the halves of a
        broken crossing (the unwound stretch has part None); transitions = the
        evidence for every fragment junction inside the crossing. AMBIGUOUS
        junctions lower the crossing's confidence one level."""
        lower = {"high": "medium", "medium": "low", "low": "low"}
        for region in regions:
            pos = frame.position(region.start, region.start_icode or "")
            if pos is None or groups[pos] < 0:
                continue
            g = groups[pos]
            region.crossing = g + 1
            if g in transitions:
                region.transitions = transitions[g]
                if region.confidence and any(t["classification"] == "AMBIGUOUS" for t in transitions[g]):
                    region.confidence = lower[region.confidence]
            if g in broken and region.description != L.UNWOUND:
                part_a = broken[g][0]
                region.part = "a" if pos <= part_a[1] else "b"

    @staticmethod
    def _attach_topology_labels(regions) -> None:
        """Membrane role + topology label for every region, e.g.
        TM1a / TM1 unwound / TM1b (BROKEN_TM, UNWOUND), TM2 (TM_CROSSING),
        EL2, EL3a / EL3b (loop helices, EXTRAMEMBRANE or INTERFACIAL), IL1,
        N-term / C-term, RE1 (REENTRANT). Secondary structure stays in `ss`."""
        n_cross = max((r.crossing or 0) for r in regions) if regions else 0
        loop_of, seen = [], 0
        for r in regions:
            if r.crossing:
                seen = max(seen, r.crossing)
                loop_of.append(None)
            else:
                loop_of.append(seen)
        names, el, il = {}, 0, 0
        for loop in sorted({x for x in loop_of if x is not None}):
            members = [r for r, x in zip(regions, loop_of) if x == loop]
            if loop == 0:
                names[loop] = ("N-term", "Nt")
            elif loop == n_cross:
                names[loop] = ("C-term", "Ct")
            else:
                side = next((r.side for r in members if r.side in L.KNOWN_SIDES), None)
                if side == L.EXTRA:
                    el += 1
                    names[loop] = (f"EL{el}",) * 2
                elif side == L.CYTO:
                    il += 1
                    names[loop] = (f"IL{il}",) * 2
                else:
                    names[loop] = (f"L{loop}",) * 2
        re_count = 0
        for loop, (label, short) in names.items():
            members = [r for r, x in zip(regions, loop_of) if x == loop]
            ss_elems = [id(r) for r in members if r.type == "Topological domain" and r.ss in ("Helix", "Strand")]
            letters = "abcdefghijklmnopqrstuvwxyz"
            for r in members:
                if r.type == L.INTRA:
                    re_count += 1
                    r.membrane_role, r.topology_label = "REENTRANT", f"RE{re_count}"
                elif r.type == L.SIGNAL:
                    r.membrane_role, r.topology_label = "SIGNAL", "Signal"
                else:
                    r.membrane_role = "INTERFACIAL" if L.INTERFACIAL in r.description else "EXTRAMEMBRANE"
                    if id(r) in ss_elems:
                        k = ss_elems.index(id(r))
                        r.topology_label = short + (letters[k] if len(ss_elems) > 1 else "")
                    else:
                        r.topology_label = label
        for r in regions:
            if not r.crossing:
                continue
            r.parent_tm = f"TM{r.crossing}"
            if r.description == L.UNWOUND:
                r.membrane_role, r.topology_label = "UNWOUND", f"TM{r.crossing} unwound"
            elif r.part:
                r.membrane_role, r.topology_label = "BROKEN_TM", f"TM{r.crossing}{r.part}"
            else:
                r.membrane_role, r.topology_label = "TM_CROSSING", f"TM{r.crossing}"

    def _residue_matrix(self, frame, base, raw_ss, coarse, membrane, labels, scores):
        rows = []
        for k, res in enumerate(frame.residues):
            depth = membrane.depth[k] if membrane is not None else None
            rows.append(ResidueAnnotation(
                index=k, residue_number=res.resseq, insertion_code=res.icode or None, aa=res.one,
                ss_raw=raw_ss[k] if raw_ss is not None else None,
                ss=coarse[k] if coarse is not None else None,
                tm_evidence=base[k],
                depth=round(float(depth), 2) if depth is not None and np.isfinite(depth) else None,
                zone=membrane.zone(k) if membrane is not None else None,
                plddt=(round(float(frame.b_factors[k]), 1)
                       if frame.predicted_model and frame.b_factors is not None else None),
                label=labels[k],
                confidence=self._level(scores[k]) if scores[k] is not None else None,
            ))
        return rows

    @staticmethod
    def _domain_type(spans) -> Optional[str]:
        classes = [c for _, _, c in spans]
        h, e = classes.count("H"), classes.count("E")
        if not spans:
            return None
        if h and not e:
            return "alpha_helical"
        if e and not h:
            return "beta_barrel" if e >= 8 else "beta"
        if h and e:
            return "mixed"
        return "irregular"

    def _validate(self, labels, spans, frame) -> list[str]:
        """Topology checks. They only report - they never change the annotation."""
        notes = []
        for i, (s0, e0, cls) in enumerate(spans, start=1):
            rng = TM_LENGTH_RANGE.get(cls)
            length = e0 - s0 + 1
            if rng and not rng[0] <= length <= rng[1]:
                kind = "helix" if cls == "H" else "strand"
                notes.append(f"validation: TM{i} ({frame.residues[s0].label}-{frame.residues[e0].label}) "
                             f"is {length} residues, outside the usual {rng[0]}-{rng[1]} for a TM {kind}")
        # sides must alternate across every crossing
        sides = []
        for (s0, e0, _), nxt in zip(spans, spans[1:] + [None]):
            stop = nxt[0] if nxt else len(labels)
            loop = [lab.split()[0] for lab in labels[e0 + 1:stop] if lab.split()[0] in L.KNOWN_SIDES]
            sides.append(loop[0] if loop else None)
        before = [lab.split()[0] for lab in labels[:spans[0][0]] if lab.split()[0] in L.KNOWN_SIDES]
        seq = [before[-1] if before else None] + sides
        for i in range(len(seq) - 1):
            if seq[i] and seq[i + 1] and seq[i] == seq[i + 1]:
                notes.append(f"validation: both sides of TM{i + 1} are {seq[i]} "
                             "(re-entrant segment, or a crossing missed next to it)")
        # positive-inside agreement on membrane-flanking loop residues
        names = frame.names
        count = {L.CYTO: [0, 0], L.EXTRA: [0, 0]}
        for s0, e0, _ in spans:
            for k in list(range(max(0, s0 - POSITIVE_INSIDE_FLANK), s0)) + \
                    list(range(e0 + 1, min(len(labels), e0 + 1 + POSITIVE_INSIDE_FLANK))):
                side = labels[k].split()[0]
                if side in count:
                    count[side][0] += names[k] in POSITIVE_RESIDUES
                    count[side][1] += 1
        dens = {side: (c[0] / c[1] if c[1] else None) for side, c in count.items()}
        d_cyto, d_extra = dens[L.CYTO], dens[L.EXTRA]
        if d_cyto is not None and d_extra is not None and d_extra > d_cyto:
            notes.append("validation: membrane-flanking Lys/Arg are denser on the Extracellular "
                         f"side ({dens[L.EXTRA]:.2f} vs {dens[L.CYTO]:.2f}) - orientation "
                         "disagrees with the positive-inside rule")
        return notes

    def _response(self, frame, params, tm_pred, name, labeler, regions, warns, membrane=None,
                  domain_type=None, residues=None) -> TopologyResponse:
        normal = membrane.normal if membrane is not None else tm_pred.membrane_normal
        return TopologyResponse(
            uniprot_id="CALCULATED",
            protein_name=name,
            gene_name="",
            organism="Computed",
            membrane_score=round(float(tm_pred.membrane_score or 0.0), 3),
            membrane_normal=[round(float(x), 4) for x in normal] if normal is not None else None,
            labeler=labeler,
            parameters_used={**params.to_response_dict(), "flow": self.flow_type,
                             "max_snap": self.max_snap},
            regions=regions,
            chain_id=frame.chain_id,
            warnings=list(dict.fromkeys(warns)),
            domain_type=domain_type,
            membrane=membrane.as_dict() if membrane is not None else None,
            residues=residues,
        )

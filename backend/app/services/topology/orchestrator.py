"""
Topology orchestrator - ONE flow (consensus).

The TM block and the SS block run side by side (concurrently) on the same residue
frame and are merged residue by residue. There are no other flows any more: the old
tm_then_ss / ss_then_tm values are accepted for API compatibility, reported in the
warnings, and run as consensus.

INDEXING CONTRACT
  The structure is parsed ONCE into a ResidueFrame (residues.py). The TM provider
  returns one label per frame POSITION, the SS provider returns one code per frame
  POSITION (mapped by residue key, verified by residue identity, alignment fallback),
  and everything below works on those position arrays only. Author residue numbers
  come back only when regions are emitted, so the two blocks can never disagree about
  which residue is which.

THE FLOW
  1. Parse the chain once -> ResidueFrame, positions 0..N-1.
  2. TM block (UniProt / Kyte-Doolittle / 3D slab) and SS block (DSSP / STRIDE) run
     at the same time; each reports one value per position.
  3. Membrane geometry (membrane.py): file planes (OPM/PPM DUM atoms) > axes of the
     TM segments. Used for depth / zone, the hairpin guard and interfacial helices -
     it never moves a TM boundary.
  4. Optional treat_turn_as_helix (see _promote_turns).
  5. Residue map (consensus.build_consensus_map), loop i = 0..N-1:
        flag_tm = residue i lies in a TM segment of the TM block (the bilayer band)
        flag_ss = residue i is helix (H/G/I; E inside a beta TM segment)
        flag_ss and flag_tm      -> TM_in
        flag_ss and not flag_tm  -> TM_C / TM_E (side of the residue)
        not flag_ss              -> not in the map
  6. Crossings are drawn from TM_in only, per TM segment of the TM block:
        - a TM_in run of ONE residue is not TM (relabelled TM_C / TM_E in the map)
        - 1-2 non-helical residues inside a run = kink (joined), unless the chain
          turns back (antiparallel along the normal = hairpin)
        - fragments shorter than MIN_FRAGMENT_LEN are not drawn as their own part
        - >= 3 non-helical residues between two runs of one segment:
             antiparallel, or each run crosses the bilayer alone -> two crossings
             otherwise -> one crossing drawn as TMa / TM unwound / TMb
        - a chain break is never called "unwound"
        - a TM segment without TM_in is not drawn (warning)
  7. Sides (TM block first, else alternation + positive-inside), loop SS, interfacial
     helices, confidence, regions, topology labels and a validation that only warns.
  Without an SS result the TM block's segments are reported unchanged (warning).
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from functools import partial
from pathlib import Path
from typing import Optional, Sequence, cast

import numpy as np

from app.core.constants import POSITIVE_RESIDUES
from app.schemas.topology import (
    ConsensusResidue, ResidueAnnotation, TMParams, TopologyRegion, TopologyResponse,
)
from app.services.topology import labels as L
from app.services.topology.consensus import (
    HELIX_CODES, TM_C, TM_E, TM_IN, TURN_IN, ConsensusEntry, build_consensus_map,
    drop_short_fragments, join_kinks, tm_in_runs,
)
from app.services.topology.membrane import MembraneFrame, build_membrane
from app.services.topology.providers.ss.base import SSProvider
from app.services.topology.providers.tm.base import TMPrediction, TMProvider
from app.services.topology.residues import (
    ResidueFrame, SSAssignment, SSRecord, load_residue_frame, map_ss_records, parse_resnum,
)
from app.services.topology.transitions import (
    MIN_FRAGMENT_LEN, TransitionEvidence, classify_membrane_transition,
)

# The only flow. Old flow names are still accepted (and reported) so existing clients
# keep working; any other value is an error.
FLOW = "consensus"
CONSENSUS_NAMES = {"consensus", "parallel_merge"}
REMOVED_FLOWS = {"tm_then_ss", "tm_first", "ss_then_tm", "structure_guided"}

# Loop residues within this distance of a TM end vote in the positive-inside rule.
POSITIVE_INSIDE_FLANK = 15
MIN_SS_COVERAGE_WARN = 0.8

Span = tuple[int, int, Optional[str]]      # (start_pos, end_pos, 'H' | 'E' | None)

# Plausible lengths of one crossing (residues) - validation warnings only.
TM_LENGTH_RANGE = {"H": (14, 40), "E": (5, 16)}
INTERFACIAL_MIN_LEN = 6
# hairpin guard: each fragment must travel at least this far along the membrane
# normal (A) for its direction to count (a 5-residue helix ~ 7.5 A along its axis)
HAIRPIN_MIN_TRAVEL = 3.0
# DSSP/STRIDE turn and bend codes (treat_turn_as_helix)
TURN_CODES = frozenset("TS")


class TopologyOrchestrator:
    def __init__(self, tm_provider: TMProvider, ss_provider: Optional[SSProvider],
                 flow_type: Optional[str] = None, *, annotate_extramembrane_ss: bool = True,
                 include_residues: bool = True):
        name = (flow_type or FLOW).strip()
        if name not in CONSENSUS_NAMES | REMOVED_FLOWS:
            raise ValueError(f"Unknown flow_type {flow_type!r}; the only flow is {FLOW!r}")
        self.flow_type = FLOW
        self._flow_notes = ([] if name in CONSENSUS_NAMES else
                            [f"flow_type '{name}' has been removed; the consensus flow was used"])
        self.tm_provider = tm_provider
        self.ss_provider = ss_provider
        self.annotate_extramembrane_ss = annotate_extramembrane_ss
        self.include_residues = include_residues

    # ------------------------------------------------------------------ public
    def execute(self, file_path: Path, params: Optional[TMParams] = None,
                chain_id: Optional[str] = None, **kwargs) -> TopologyResponse:
        params = params if params is not None else TMParams()
        file_path = Path(file_path)
        frame = load_residue_frame(file_path, chain_id)
        warns: list[str] = list(self._flow_notes)
        if len(frame) == 0:
            return self._response(frame, params, TMPrediction(), "No amino-acid residues",
                                  self._tm_name(), [], warns + ["no amino-acid residues in the chain"])

        # 1. TM block and SS block, side by side on the same frame
        tm_pred, (coarse, raw_ss, ss_name, ss_warns) = self._run_blocks(file_path, frame, params,
                                                                        kwargs)
        warns.extend(tm_pred.warnings)
        base, segments, tm_warns = self._tm_on_frame(tm_pred, frame)
        warns.extend(tm_warns)
        warns.extend(ss_warns)
        if not segments:
            return self._response(frame, params, tm_pred,
                                  f"No transmembrane segments ({tm_pred.labeler or self._tm_name()})",
                                  tm_pred.labeler or self._tm_name(), [], warns)

        # 2. membrane geometry (depth / zone / hairpin guard - never boundaries)
        breaks = frame.chain_breaks()
        membrane = build_membrane(frame, segments, params.membrane_thickness / 2.0)
        if membrane is None:
            warns.append("no membrane geometry (coordinates missing): hairpins inside one TM "
                         "segment cannot be detected")
        elif coarse is not None:
            # centre an estimated bilayer on the long helices/strands of the TM segments
            spanning = [(a, b) for a, b, _ in self._ss_elements(coarse, breaks)
                        if b - a + 1 >= 8 and any(a <= e and b >= s for s, e in segments)
                        and membrane.span(a, b) >= params.min_cross_span_frac * 2 * membrane.half_thickness]
            membrane = membrane.recentered(spanning)
        labeler = f"{tm_pred.labeler or self._tm_name()} + {ss_name} ({FLOW})"

        # 3. residue map -> crossings
        consensus = None                    # (codes used by the map, strand segments, turns)
        broken: dict = {}
        transitions: dict = {}
        if coarse is None or raw_ss is None:
            spans: list[Span] = [(s, e, None) for s, e in segments]
            warns.append("no SS result: the TM block's segments are reported without the "
                         "residue-level consensus")
        else:
            map_codes, turn_indices = list(raw_ss), frozenset()
            if getattr(params, "treat_turn_as_helix", False):
                map_codes, ss_codes, turn_indices = self._promote_turns(raw_ss, breaks, membrane)
                coarse = [L.coarse_ss(c) for c in ss_codes]
            spans, broken, transitions, strand, notes = self._flow_consensus(
                base, segments, map_codes, coarse, breaks, frame, membrane, params, ss_name,
                turn_indices)
            warns.extend(notes)
            consensus = (map_codes, strand, turn_indices)
        if not spans:
            return self._response(frame, params, tm_pred,
                                  "No transmembrane segment supported by the evidence",
                                  labeler, [], warns, membrane=membrane)

        # 4. paint, sides, loop annotation
        for idx, evs in transitions.items():
            if any(ev["classification"] == "AMBIGUOUS" for ev in evs):
                warns.append(f"TM{idx + 1}: one fragment crosses the bilayer alone and the next "
                             "only partly - drawn as one crossing, confidence lowered")
        labels, groups = self._paint(base, spans)
        for _, (u0, u1), _ in broken.values():
            for k in range(u0, u1 + 1):
                labels[k] = L.UNWOUND
        labels = self._resolve_sides(labels, spans, frame)
        if coarse is not None and membrane is not None:
            labels = self._mark_interfacial(labels, coarse, breaks, frame.coords, membrane)
        if coarse is not None and self.annotate_extramembrane_ss:
            breaks_list = cast(Sequence[bool], breaks.tolist())
            labels = L.label_extramembrane_ss(labels, coarse, breaks=breaks_list)

        # 5. evidence matrix, confidence, regions, validation, residue map
        conf = self._confidence(labels, spans, base, coarse, membrane, frame)
        regions = L.build_regions(frame, labels, groups)
        self._attach_region_confidence(regions, frame, conf)
        self._attach_crossings(regions, frame, groups, broken, transitions)
        self._attach_topology_labels(regions)
        warns.extend(self._validate(labels, spans, frame))
        residues = (self._residue_matrix(frame, base, raw_ss, coarse, membrane, labels, conf)
                    if self.include_residues else None)
        consensus_map = None
        if consensus is not None:
            map_codes, strand, turn_indices = consensus
            consensus_map, notes = self._consensus_map_rows(
                frame, base, map_codes, raw_ss, labels, groups, segments, strand, turn_indices,
                broken, breaks)
            warns.extend(notes)
        return self._response(frame, params, tm_pred, "Predicted Topology", labeler, regions,
                              warns, membrane=membrane, domain_type=self._domain_type(spans),
                              residues=residues, consensus_map=consensus_map)

    # ------------------------------------------------------- blocks -> frame
    def _run_blocks(self, file_path: Path, frame: ResidueFrame, params: TMParams, kwargs):
        """TM block and SS block at the same time (UniProt download and the DSSP/STRIDE
        binary are both I/O bound). Both only READ the shared frame. An error of the TM
        block propagates (HTTP 4xx/5xx from the provider); an SS failure is a warning."""
        with ThreadPoolExecutor(max_workers=2) as pool:
            tm_future = pool.submit(self.tm_provider.predict_tm, file_path, params=params,
                                    frame=frame, **kwargs)
            ss_future = pool.submit(self._ss_on_frame, file_path, frame)
            return tm_future.result(), ss_future.result()

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
        raw_segments = tm_pred.segments if use_given and tm_pred.segments is not None else painted_segments
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

    # --------------------------------------------------------- the one flow
    @staticmethod
    def _majority(codes: Sequence[Optional[str]]) -> str:
        """Dominant regular class; 'I' when there is no helix/strand at all."""
        h, e = codes.count("H"), codes.count("E")
        if h == 0 and e == 0:
            return "I"
        return "H" if h >= e else "E"

    @staticmethod
    def _antiparallel(membrane: MembraneFrame, r1: tuple[int, int], r2: tuple[int, int]) -> bool:
        """The chain runs back across the membrane: r1 and r2 travel in opposite
        directions along the normal (both long enough to have a direction)."""
        if min(r1[1] - r1[0], r2[1] - r2[0]) + 1 < MIN_FRAGMENT_LEN:
            return False
        d = membrane.depth
        du, dv = d[r1[1]] - d[r1[0]], d[r2[1]] - d[r2[0]]
        if not (np.isfinite(du) and np.isfinite(dv)):
            return False
        return bool(du * dv < 0 and min(abs(du), abs(dv)) >= HAIRPIN_MIN_TRAVEL)

    def _promote_turns(self, raw: Sequence[Optional[str]], breaks,
                       membrane: Optional[MembraneFrame]):
        """treat_turn_as_helix. For every run of T/S residues (never across a chain break):
             helix on BOTH sides  -> the run is helix (H) - a kink inside one helix;
                                     but if the two flanking helices run antiparallel
                                     across the membrane it is a real hairpin turn and
                                     is handled as below instead
             helix on ONE side    -> Turn (Turn_in / Turn_C / Turn_E in the map)
             no helix beside it   -> unchanged
        Returns (codes for the residue map, codes for SS everywhere else, turn positions).
        Turn residues are 'H' only for the map (so they get a Turn_* entry); elsewhere
        they keep T/S, so loop helices and crossings are not stretched by turns."""
        n = len(raw)
        for_map, for_ss = list(raw), list(raw)
        turns: set[int] = set()
        i = 0
        while i < n:
            if raw[i] not in TURN_CODES:
                i += 1
                continue
            j = i + 1
            while j < n and raw[j] in TURN_CODES and not breaks[j]:
                j += 1
            left = right = None
            if i > 0 and raw[i - 1] in HELIX_CODES and not breaks[i]:
                a = i - 1
                while a > 0 and raw[a - 1] in HELIX_CODES and not breaks[a]:
                    a -= 1
                left = (a, i - 1)
            if j < n and raw[j] in HELIX_CODES and not breaks[j]:
                b = j
                while b + 1 < n and raw[b + 1] in HELIX_CODES and not breaks[b + 1]:
                    b += 1
                right = (j, b)
            hairpin = (left is not None and right is not None and membrane is not None
                       and self._antiparallel(membrane, left, right))
            if left is not None and right is not None and not hairpin:
                for k in range(i, j):
                    for_map[k] = for_ss[k] = "H"
            elif left is not None or right is not None:
                for k in range(i, j):
                    for_map[k] = "H"
                    turns.add(k)
            i = j
        return for_map, for_ss, frozenset(turns)

    def _flow_consensus(self, base, segments, map_codes, coarse, breaks, frame: ResidueFrame,
                        membrane: Optional[MembraneFrame], params: TMParams, ss_name: str,
                        turn_indices: frozenset[int] = frozenset()):
        """Residue map TM block x SS block (consensus.build_consensus_map) -> crossings.
        Crossings are drawn from TM_in residues only.
        Returns (spans, broken parts, transitions, strand segments, notes)."""
        n = len(base)
        notes: list[str] = []
        # beta-barrel TM segments: flag_ss must look for strands there, not helices
        strand = frozenset(k for k, (s, e) in enumerate(segments)
                           if self._majority(coarse[s:e + 1]) == "E")
        cmap, _ = build_consensus_map(base, map_codes, tm_segments=segments,
                                      strand_segments=strand, turn_indices=turn_indices)
        # a TM_in run of a single residue is not TM (it becomes TM_C / TM_E in the map)
        runs = {k: [r for r in rs if r[1] > r[0]] for k, rs in tm_in_runs(cmap, n, breaks).items()}
        is_hairpin = partial(self._antiparallel, membrane) if membrane is not None else None
        min_cross = full_cross = 0.0
        if membrane is not None:
            thickness = 2.0 * membrane.half_thickness
            min_cross = params.min_cross_span_frac * thickness
            full_cross = params.full_cross_frac * thickness

        # crossings: (members, junctions, class); junction = ("unwound", evidence | None)
        # or ("chain_break", None)
        crossings = []
        for k, (s0, e0) in enumerate(segments):
            cls = "E" if k in strand else "H"
            tm_name = f"{frame.residues[s0].label}-{frame.residues[e0].label}"
            seg_runs = drop_short_fragments(join_kinks(runs.get(k, []), breaks,
                                                       is_hairpin=is_hairpin))
            if not seg_runs:
                notes.append(f"TM segment {tm_name} of the TM block has no "
                             f"{'strand' if cls == 'E' else 'helix'} residue in {ss_name} "
                             "(no TM_in) - not drawn")
                continue
            members = [seg_runs[0]]
            junctions: list[tuple[str, Optional[TransitionEvidence]]] = []
            for run in seg_runs[1:]:
                prev = members[-1]
                if any(breaks[p] for p in range(prev[1] + 1, run[0] + 1)):
                    # residues missing from the model: unresolved is not "unwound"
                    members.append(run)
                    junctions.append(("chain_break", None))
                    notes.append(f"chain break inside TM segment {tm_name} - drawn as one "
                                 "crossing, not evaluated as a broken helix")
                    continue
                ev = None
                if membrane is not None:
                    ev = classify_membrane_transition((members[0][0], prev[1]), run, membrane,
                                                      frame.coords, breaks, min_cross, full_cross)
                    antiparallel = self._antiparallel(membrane, prev, run)
                    if ev.classification == "TWO_TM" or antiparallel:
                        # two helices under one TM feature (typical of sequence-predicted
                        # or slab-fit TM segments): two crossings, never TMa/TMb
                        crossings.append((members, junctions, cls))
                        members, junctions = [run], []
                        notes.append(f"TM segment {tm_name} holds two "
                                     f"{'antiparallel ' if antiparallel else ''}helices "
                                     f"({ev.reason}) - drawn as two crossings")
                        continue
                members.append(run)
                junctions.append(("unwound", ev))
            crossings.append((members, junctions, cls))

        crossings.sort(key=lambda c: c[0][0][0])
        spans: list[Span] = []
        broken: dict = {}
        transitions: dict = {}
        for idx, (members, junctions, cls) in enumerate(crossings):
            s, e = members[0][0], members[-1][1]
            spans.append((s, e, cls))
            evs: list[dict] = []
            for kind, ev in junctions:
                if ev is not None:
                    d = ev.as_dict()
                    d["decision"] = "BROKEN_TM"
                    d["decision_source"] = "consensus map: non-helical residues inside the TM band"
                    evs.append(d)
            if evs:
                transitions[idx] = evs
            cuts = [(members[i + 1][0] - members[i][1] - 1, i)
                    for i, (kind, _) in enumerate(junctions) if kind == "unwound"]
            if not cuts:
                continue
            _, i = max(cuts)
            left_end, right_start = members[i][1], members[i + 1][0]
            broken[idx] = ((s, left_end), (left_end + 1, right_start - 1), (right_start, e))
            if len(cuts) > 1:
                notes.append(f"TM{idx + 1} has {len(cuts)} unwound stretches; only the widest "
                             "is drawn as the a/b break")
        return spans, broken, transitions, strand, notes

    @staticmethod
    def _nearest_side(labels: Sequence[str], i: int) -> Optional[str]:
        """Side (Cytoplasmic / Extracellular) of the nearest residue that has one;
        on a tie the N-terminal neighbour wins. Labels may carry an SS suffix
        ("Cytoplasmic Helix"), so only the first word is read."""
        n = len(labels)
        for step in range(1, n):
            for k in (i - step, i + step):
                if 0 <= k < n and labels[k].split()[0] in L.KNOWN_SIDES:
                    return labels[k].split()[0]
        return None

    def _consensus_map_rows(self, frame: ResidueFrame, base, map_codes, raw_ss, labels, groups,
                            segments, strand, turn_indices: frozenset[int], broken: dict, breaks):
        """Final residue map with the resolved sides (TM block's side first, inferred
        otherwise) -> response rows + notes."""
        cmap, skipped = build_consensus_map(base, map_codes, sides=labels, tm_segments=segments,
                                            strand_segments=strand, turn_indices=turn_indices)
        # a TM_in run of ONE residue is not TM: it takes the side of the nearest residue
        for rs in tm_in_runs(cmap, len(frame), breaks).values():
            for s, e in rs:
                if s == e:
                    side = self._nearest_side(labels, s)
                    if side is not None:
                        cmap[s] = ConsensusEntry(TM_C if side == L.CYTO else TM_E,
                                                 cmap[s].ss_raw, None)
        rows = []
        for i in sorted(cmap):
            entry, res = cmap[i], frame.residues[i]
            g = groups[i]
            crossing = g + 1 if entry.label in (TM_IN, TURN_IN) and g >= 0 else None
            part = None
            if crossing is not None and g in broken:
                (_, a_end), _, (b_start, _) = broken[g]
                part = "a" if i <= a_end else "b" if i >= b_start else None
            rows.append(ConsensusResidue(
                index=i, residue_number=res.resseq, insertion_code=res.icode or None, aa=res.one,
                label=entry.label,
                ss_raw=raw_ss[i] if raw_ss[i] is not None else entry.ss_raw,
                tm_segment=entry.tm_segment + 1 if entry.tm_segment is not None else None,
                crossing=crossing, part=part,
            ))
        notes = []
        if skipped:
            what = ", ".join(f"{v} in {k}" for k, v in sorted(skipped.items()))
            notes.append(f"consensus map: helix residues outside the TM band without a "
                         f"cytoplasmic/extracellular side get no TM_C/TM_E label ({what})")
        return rows, notes

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
            idx: set[int] = set()
            if a > 0:                                   # TM span before the loop
                idx.update(range(a, min(b, a + POSITIVE_INSIDE_FLANK - 1) + 1))
            if b < n - 1:                               # TM span after the loop
                idx.update(range(max(a, b - POSITIVE_INSIDE_FLANK + 1), b + 1))
            pos[parity] += sum(names[p] in POSITIVE_RESIDUES for p in idx)
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
        loop_of: list[Optional[int]] = []
        seen = 0
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
        dens_cyto = dens.get(L.CYTO)
        dens_extra = dens.get(L.EXTRA)
        if dens_cyto is not None and dens_extra is not None and dens_extra > dens_cyto:
            notes.append("validation: membrane-flanking Lys/Arg are denser on the Extracellular "
                         f"side ({dens_extra:.2f} vs {dens_cyto:.2f}) - orientation "
                         "disagrees with the positive-inside rule")
        return notes

    def _response(self, frame, params, tm_pred, name, labeler, regions, warns, membrane=None,
                  domain_type=None, residues=None, consensus_map=None) -> TopologyResponse:
        normal = membrane.normal if membrane is not None else tm_pred.membrane_normal
        return TopologyResponse(
            uniprot_id="CALCULATED",
            protein_name=name,
            gene_name="",
            organism="Computed",
            membrane_score=round(float(tm_pred.membrane_score or 0.0), 3),
            membrane_normal=[round(float(x), 4) for x in normal] if normal is not None else None,
            labeler=labeler,
            parameters_used={**params.to_response_dict(), "flow": FLOW},
            regions=regions,
            chain_id=frame.chain_id,
            warnings=list(dict.fromkeys(warns)),
            domain_type=domain_type,
            membrane=membrane.as_dict() if membrane is not None else None,
            residues=residues,
            consensus_map=consensus_map,
        )

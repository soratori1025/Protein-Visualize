"""
"3D Energy" TM provider: the membrane is placed by minimising an implicit-membrane
transfer energy (see `app.services.topology.membrane_energy`), then the analysed chain is
labelled from that placement.

Differences from the 3D slab (`GeometryTMProvider`, which is kept unchanged):
  * only solvent-exposed surface counts (Shrake-Rupley SASA of the whole model), so the
    buried hydrophobic core of a soluble domain and subunit interfaces play no role;
  * residues lining a water-filled pore are discounted (ray-escape test);
  * a depth profile with a hydrocarbon core and an interface region (Trp/Tyr belt)
    replaces the hard 30 A slab;
  * the hydrophobic thickness is fitted (20-40 A) together with the normal and centre;
  * every protein chain of the model takes part in the fit; the chain asked for is
    labelled;
  * the membrane/soluble decision is a transfer free energy in kcal/mol;
  * predicted models: residues with pLDDT < 50 are left out of the fit, and in-membrane
    runs made mostly of them are not reported as TM;
  * proteins that bend the bilayer (Piezo) get a spherical membrane (PPM 3.0 style) when
    it fits clearly better than a plane; the fitted geometry (depth of every residue,
    thickness, curvature) is returned so the orchestrator does not re-estimate a plane.

Labels use the same post-processing as the slab (chain breaks, jitter vs. turns, split at
turns, short / non-crossing / hydrophilic runs) so both providers are comparable.
`membrane_score` holds -dG_transfer (kcal/mol, larger = more membrane-like).
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np

from app.schemas.topology import TMParams
from app.services.topology import labels as L
from app.services.topology import membrane_energy as ME
from app.services.topology.providers.tm.base import TMPrediction, TMProvider
from app.services.topology.residues import ResidueFrame, load_residue_frame


def _drop_low_confidence(classes, frame, d, T, predicted: bool = False):
    """Predicted models: an in-membrane run made mostly of residues with pLDDT <
    MIN_PLDDT is a low-confidence loop that happens to cross the bilayer, not a TM."""
    b = getattr(frame, "b_factors", None)
    predicted = predicted or bool(getattr(frame, "predicted_model", False))
    if not predicted or b is None or len(b) != len(classes):
        return classes, 0
    b = np.asarray(b, dtype=float)
    if np.nanmax(b) <= 1.0:
        b = b * 100.0
    out, dropped = list(classes), 0
    for s, e in T._runs_where([c == T.TM_LABEL for c in out], T._chain_breaks(frame.coords)):
        if np.mean(b[s:e + 1] < ME.MIN_PLDDT) > 0.5:
            side = T._side_for(np.asarray(d)[s:e + 1])
            for k in range(s, e + 1):
                out[k] = side
            dropped += 1
    return out, dropped


class EnergyTMProvider(TMProvider):
    LABELER = "3D_Energy"

    @staticmethod
    def _placement_quality(half_thickness: float, energy: float, n_segments: int,
                            n_chains: int, predicted: bool):
        """Is the placement determined by the data, or is it a guess?

        Three situations make the fit under-determined, and all three have been
        observed:

        * the fitted hydrophobic thickness sits at the edge of the search range,
          which happens when the structure does not constrain the bilayer - a
          monomeric subunit of a curved assembly, for instance;
        * the transfer free energy is barely below the decision threshold, so the
          structure is only marginally membrane-like;
        * very few transmembrane segments were found in a large assembly, so the
          membrane is inferred from a handful of helices.

        A predicted structure is reported separately: the fact that the model is a
        prediction is a caveat, not a defect of the placement, and it must not
        drown out the geometric verdict.

        The result is a set of flags, a set of notes and a coarse quality label.
        """
        flags = []
        notes = []
        tolerance = 0.05 * ME.MAX_HALF
        if half_thickness <= ME.MIN_HALF + tolerance:
            flags.append("the fitted hydrophobic thickness is at the lower bound of the "
                         "search range: the bilayer is not well determined by this structure")
        elif half_thickness >= ME.MAX_HALF - tolerance:
            flags.append("the fitted hydrophobic thickness is at the upper bound of the "
                         "search range: the bilayer is not well determined by this structure")
        if energy > 1.5 * ME.MAX_MEMBRANE_DG:
            flags.append(f"the transfer free energy ({energy:.1f} kcal/mol) is only just "
                         f"below the membrane threshold: the structure is marginally "
                         f"membrane-like")
        if n_segments <= 2 and n_chains > 1:
            flags.append(f"only {n_segments} transmembrane segment(s) were found in an "
                         f"assembly of {n_chains} chains: the membrane placement rests on "
                         f"very little evidence")
        if predicted:
            notes.append("the structure is a prediction: the placement depends on the "
                         "pLDDT filter, which excludes the low-confidence residues")
        if not flags:
            return flags, notes, "ok"
        if any("not well determined" in f or "very little evidence" in f for f in flags):
            return flags, notes, "uncertain"
        return flags, notes, "caution"

    def predict_tm(self, file_path: Path, params: Optional[TMParams] = None,
                   frame: Optional[ResidueFrame] = None, **kwargs) -> TMPrediction:
        # shared post-processing lives with the slab provider (imported lazily: large module)
        from app.services.topology import topology_predictor as T

        frame = frame if frame is not None else load_residue_frame(file_path, kwargs.get("chain_id"))
        n = len(frame)
        if n < 10:
            return TMPrediction(labeler=self.LABELER,
                                warnings=[f"only {n} residues - too few to place a membrane"])
        coords = np.asarray(frame.coords, dtype=float)
        fit = ME.fit_structure(file_path, ca_coords=coords, ca_names=frame.names)
        pl = fit.placement
        warnings = list(fit.warnings)
        dg, half = pl.energy, pl.half_thickness
        labeler = (
            f"{self.LABELER} (dG {dg:.1f} kcal/mol, hydrophobic thickness "
            f"{2 * half:.1f} A, {fit.mode}"
            + (f", curved R={pl.radius:.0f} A" if (pl.curved and pl.radius is not None) else "")
            + (f", {fit.n_chains} chains" if fit.n_chains > 1 else "")
            + ")"
        )
        normal = [round(float(a), 4) for a in pl.normal]
        d = pl.depth(coords)
        residues_data = frame.residues_data()
        geometry = dict(  # handed to the orchestrator's membrane layer
            depth=[round(float(x), 3) for x in d],
            half_thickness=float(half),
            membrane_center=[float(x) for x in pl.center],
            sphere_center=[float(x) for x in pl.sphere_center] if (pl.curved and pl.sphere_center is not None) else None,
            membrane_radius=float(pl.radius) if (pl.curved and pl.radius is not None) else None,
        )

        if dg > ME.MAX_MEMBRANE_DG:
            labels = T._apply_positive_inside_rule(T._sides_by_sign(d), residues_data)
            return TMPrediction(
                labels=labels, segments=[], membrane_score=-dg, membrane_normal=normal,
                labeler=labeler, **geometry,
                placement_quality="uncertain",
                placement_flags=[f"the transfer free energy ({dg:.1f} kcal/mol) is above "
                                 f"the membrane threshold, so the structure is treated as "
                                 f"soluble"],
                warnings=warnings + [f"dG_transfer {dg:.1f} kcal/mol > {ME.MAX_MEMBRANE_DG} "
                                     "- treated as soluble (no TM segments)"])

        breaks = T._chain_breaks(coords)
        raw_kd, _ = T._kd_raw(frame.names)
        classes = [T.TM_LABEL if abs(x) <= half else ("Side_A" if x > 0 else "Side_B") for x in d]
        classes = T._smooth_flickers(classes, d, half, breaks=breaks)
        classes = T._split_at_turns(classes, d, half, breaks=breaks)
        classes = T._drop_short_tm(classes, d=d, breaks=breaks)
        classes, n_dropped = T._drop_non_crossing(classes, d, half, raw_kd, breaks=breaks)
        classes, n_low = _drop_low_confidence(classes, frame, d, T, fit.predicted)
        if n_low:
            warnings.append(f"{n_low} in-membrane run(s) dropped: mostly pLDDT < "
                            f"{ME.MIN_PLDDT:.0f} in this predicted model")
        if n_dropped:
            warnings.append(f"{n_dropped} in-membrane run(s) dropped: they do not cross the "
                            "membrane or are strongly hydrophilic")
        if not any(c == T.TM_LABEL for c in classes):
            warnings.append("the membrane was placed on the assembly but this chain has no "
                            "residue crossing it")
        labels = T._apply_positive_inside_rule(classes, residues_data)
        segments = T._tm_runs(labels, breaks=breaks)
        flags, notes, quality = self._placement_quality(half, dg, len(segments),
                                                        fit.n_chains, fit.predicted)
        for flag in flags:
            warnings.append(f"placement quality: {flag}")
        for note in notes:
            warnings.append(f"placement note: {note}")
        return TMPrediction(
            regions=L.build_regions(frame, labels, T._groups_for_runs(n, segments)),
            labels=labels,
            segments=segments,
            membrane_normal=normal,
            membrane_score=-dg,
            labeler=labeler,
            warnings=warnings,
            placement_flags=flags + notes,
            placement_quality=quality,
            **geometry,
        )
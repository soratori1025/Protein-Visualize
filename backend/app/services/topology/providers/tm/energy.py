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
  * the membrane/soluble decision is a transfer free energy in kcal/mol.

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


class EnergyTMProvider(TMProvider):
    LABELER = "3D_Energy"

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
        labeler = (f"{self.LABELER} (dG {dg:.1f} kcal/mol, hydrophobic thickness "
                   f"{2 * half:.1f} A, {fit.mode}"
                   + (f", {fit.n_chains} chains" if fit.n_chains > 1 else "") + ")")
        normal = [round(float(a), 4) for a in pl.normal]
        d = pl.depth(coords)
        residues_data = frame.residues_data()

        if dg > ME.MAX_MEMBRANE_DG:
            labels = T._apply_positive_inside_rule(T._sides_by_sign(d), residues_data)
            return TMPrediction(
                labels=labels, segments=[], membrane_score=-dg, membrane_normal=normal,
                labeler=labeler,
                warnings=warnings + [f"dG_transfer {dg:.1f} kcal/mol > {ME.MAX_MEMBRANE_DG} "
                                     "- treated as soluble (no TM segments)"])

        breaks = T._chain_breaks(coords)
        raw_kd, _ = T._kd_raw(frame.names)
        classes = [T.TM_LABEL if abs(x) <= half else ("Side_A" if x > 0 else "Side_B") for x in d]
        classes = T._smooth_flickers(classes, d, half, breaks=breaks)
        classes = T._split_at_turns(classes, d, half, breaks=breaks)
        classes = T._drop_short_tm(classes, d=d, breaks=breaks)
        classes, n_dropped = T._drop_non_crossing(classes, d, half, raw_kd, breaks=breaks)
        if n_dropped:
            warnings.append(f"{n_dropped} in-membrane run(s) dropped: they do not cross the "
                            "membrane or are strongly hydrophilic")
        if not any(c == T.TM_LABEL for c in classes):
            warnings.append("the membrane was placed on the assembly but this chain has no "
                            "residue crossing it")
        labels = T._apply_positive_inside_rule(classes, residues_data)
        segments = T._tm_runs(labels, breaks=breaks)
        return TMPrediction(
            regions=L.build_regions(frame, labels, T._groups_for_runs(n, segments)),
            labels=labels,
            segments=segments,
            membrane_normal=normal,
            membrane_score=-dg,
            labeler=labeler,
            warnings=warnings,
        )
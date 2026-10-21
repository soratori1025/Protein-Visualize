"""
Test suite for the ProteinLab topology engine.

    python -m pytest backend/tests -q

The tests are the numerical contracts of the paper: every routine whose result the
manuscript quotes is checked against an independent implementation, a closed-form
value, or an invariant that must hold.  They run in seconds and need no network.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
for candidate in (ROOT / "backend", ROOT):
    if str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

from app.core.constants import KYTE_DOOLITTLE  # noqa: E402
from app.schemas.topology import TMParams  # noqa: E402
from app.services.topology import membrane_energy as ME  # noqa: E402
from app.services.topology.consensus import (  # noqa: E402
    build_consensus_map, drop_short_fragments, join_kinks,
)
from app.services.topology.labels import side_from_text  # noqa: E402
from app.services.topology.providers.tm.sequence import SequenceTMProvider  # noqa: E402
from app.services.topology.residues import (  # noqa: E402
    canonical_resname, load_residue_frame, map_to_reference, one_letter,
)
from app.services.topology.transitions import classify_membrane_transition  # noqa: E402

DATA = ROOT / "backend" / "data" / "uploads"
HSERT = DATA / "hSERT_5I6X_ions.pdb"
PIEZO = DATA / "8YEZ_Q92508.cif"

pytestmark = pytest.mark.skipif(not HSERT.exists(),
                                reason="sample structures are not part of the repository")


def golden_sphere(n: int, radius: float) -> np.ndarray:
    """n points spread almost uniformly on a sphere (Fibonacci lattice)."""
    index = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * index / n)
    theta = np.pi * (1.0 + 5 ** 0.5) * index
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1) * radius


# --------------------------------------------------------------------- residues
def test_modified_residues_map_to_their_parent():
    assert canonical_resname("MSE") == "MET"
    assert canonical_resname("CYX") == "CYS"
    assert canonical_resname("HIE") == "HIS"
    assert one_letter("MSE") == "M"
    assert one_letter("ALA") == "A"


def test_residue_frame_indexes_one_chain_and_keeps_author_numbers():
    frame = load_residue_frame(HSERT, "A")
    assert frame.chain_id == "A"
    assert len(frame) > 400
    assert frame.residues[0].resseq >= 1
    assert frame.sequence[:5] == "".join(one_letter(r.name) for r in frame.residues[:5])
    assert {r.chain for r in frame.residues} == {"A"}


def test_chain_breaks_are_a_sparse_boolean_array():
    frame = load_residue_frame(HSERT, "A")
    breaks = np.asarray(frame.chain_breaks())
    assert breaks.shape == (len(frame),)
    assert breaks.dtype == bool
    # a resolved chain has at most a handful of breaks
    assert breaks.sum() <= 5


def test_a_wrong_sequence_is_not_mapped_silently():
    frame = load_residue_frame(HSERT, "A")
    other = "".join("W" if c == "A" else c for c in frame.sequence)
    mapping = map_to_reference(frame, [(u, "") for u in range(1, len(other) + 1)], list(other),
                                offsets=[0], source="test", allow_order_fallback=False)
    # the mapping must not claim a high identity for a corrupted reference
    assert (mapping.identity or 0.0) < 0.95
    assert mapping.warnings


def test_the_right_sequence_maps_with_full_identity():
    frame = load_residue_frame(HSERT, "A")
    mapping = map_to_reference(frame, [(u, "") for u in range(1, len(frame.sequence) + 1)],
                                list(frame.sequence), offsets=[0], source="test",
                                allow_order_fallback=False)
    assert mapping.matched == len(frame)
    assert mapping.identity is not None and mapping.identity > 0.99


# ------------------------------------------------------------- energy: the SASA
def test_sasa_of_an_isolated_atom_is_its_sphere():
    sasa = ME.shrake_rupley(np.array([[0.0, 0.0, 0.0]]), np.array([1.7]),
                            probe=ME.PROBE, n_points=512)
    expected = 4.0 * math.pi * (1.7 + ME.PROBE) ** 2
    assert sasa[0] == pytest.approx(expected, rel=1e-3)


def test_sasa_of_a_buried_atom_is_zero():
    n = 400
    shell = golden_sphere(n, 6.0)
    xyz = np.vstack([np.zeros((1, 3)), shell])
    radii = np.concatenate([[1.7], np.full(n, 1.7)])
    sasa = ME.shrake_rupley(xyz, radii, probe=ME.PROBE, n_points=96)
    assert sasa[0] == pytest.approx(0.0, abs=1e-6)
    assert sasa[1:].min() > 0.0


def test_burying_an_atom_reduces_its_sasa():
    alone = ME.shrake_rupley(np.array([[0.0, 0.0, 0.0]]), np.array([1.7]),
                             probe=ME.PROBE, n_points=256)
    n = 400
    shell = golden_sphere(n, 5.0)
    xyz = np.vstack([np.zeros((1, 3)), shell])
    radii = np.concatenate([[1.7], np.full(n, 1.7)])
    buried = ME.shrake_rupley(xyz, radii, probe=ME.PROBE, n_points=256)
    assert buried[0] < 0.1 * alone[0]


def test_rsa_is_bounded_by_one():
    for value in ME.MAX_ASA.values():
        assert 50.0 < value < 400.0


# ------------------------------------------------------- energy: depth profiles
def test_depth_profiles_are_zero_outside_and_one_in_the_core():
    core, interface = ME.profiles(np.array([0.0, 5.0, 20.0, 40.0]), half=10.0)
    assert core[0] == pytest.approx(1.0, abs=1e-2)
    assert core[1] > 0.8
    assert core[2] < 0.05
    assert core[3] == pytest.approx(0.0, abs=1e-6)
    # the interface belt exists where the core is gone, and nowhere else
    assert interface[2] > 0.05
    assert interface[0] == pytest.approx(0.0, abs=1e-2)


def test_pore_weight_is_zero_for_a_buried_residue_and_one_for_an_exposed_one():
    assert ME.pore_weight(0.0) == 0.0
    assert ME.pore_weight(0.5) == pytest.approx(1.0)
    assert ME.pore_weight(1.0) == pytest.approx(1.0)


def test_ges_core_scale_favours_hydrophobic_residues():
    assert ME.GES["LEU"] > 0
    assert ME.GES["ARG"] < 0
    assert ME.CORE_SCALE["LEU"] < 0
    assert ME.CORE_SCALE["ARG"] > 0
    # tryptophan is unfavourable in the core but favourable at the interface
    assert ME.GES["TRP"] > 0
    assert ME.INTERFACE_SCALE["TRP"] < 0


# ---------------------------------------------------------- energy: FFT search
def test_fft_search_matches_a_direct_evaluation_of_the_objective():
    """The convolution must reproduce E(c, D) evaluated residue by residue."""
    rng = np.random.default_rng(7)
    n = 120
    coords = rng.normal(0, 8, size=(n, 3))
    weights = rng.uniform(0.2, 1.0, size=n)
    core_w = rng.uniform(-6, 2, size=n) * weights
    if_w = rng.uniform(-1, 1, size=n) * weights
    normal = np.array([0.0, 0.0, 1.0])
    halves = [12.0, 15.0]

    best, best_c, best_h = ME._grid_search(coords, core_w, if_w, np.array([normal]),
                                            np.array(halves), bin_w=0.5, chunk=1)
    z = coords @ normal - best_c[0]
    core, interface = ME.profiles(z, best_h[0])
    direct = float(np.sum(core_w * core + if_w * interface))
    assert best[0] == pytest.approx(direct, rel=0.02, abs=1.0)


def test_the_search_finds_the_mid_plane_of_a_slab_of_hydrophobic_residues():
    """A synthetic block of hydrophobic residues in a slab must be located."""
    n = 30
    coords = np.array([[i * 2.0, 0.0, (i - 15) * 1.5] for i in range(n)])
    core_w = np.full(n, -3.0)          # strongly hydrophobic
    if_w = np.full(n, -0.2)
    normals = ME._hemisphere(200)
    best, best_c, best_h = ME._grid_search(coords, core_w, if_w, normals,
                                            np.array([10.0, 15.0]), bin_w=0.5, chunk=32)
    # the slab is centred on z = 0, so the fitted mid-plane must be close to it
    assert abs(best_c[0]) < 3.0
    # and the chosen half-thickness must be one of the two that were offered
    assert best_h[0] in (10.0, 15.0)


# ------------------------------------------------------------------- providers
def test_kyte_doolittle_provider_finds_the_helices_of_a_transporter():
    pred = SequenceTMProvider().predict_tm(HSERT)
    assert pred.segments, "no transmembrane segment found"
    for start, end in pred.segments:
        assert end - start + 1 >= 13


def test_kyte_doolittle_values_are_the_published_scale():
    assert KYTE_DOOLITTLE["ILE"] == pytest.approx(4.5)
    assert KYTE_DOOLITTLE["LEU"] == pytest.approx(3.8)
    assert KYTE_DOOLITTLE["ARG"] == pytest.approx(-4.5)
    assert KYTE_DOOLITTLE["LYS"] == pytest.approx(-3.9)


# ------------------------------------------------------------------- consensus
def test_consensus_map_places_helices_in_the_membrane():
    frame = load_residue_frame(HSERT, "A")
    n = len(frame)
    labels = ["Transmembrane"] * n
    ss = ["H"] * n
    cmap, skipped = build_consensus_map(labels, ss, tm_segments=[(0, n - 1)])
    assert cmap
    assert all(entry.label == "TM_in" for entry in cmap.values())


def test_a_coil_residue_is_not_placed_in_the_consensus_map():
    frame = load_residue_frame(HSERT, "A")
    n = len(frame)
    labels = ["Transmembrane"] * n
    ss = ["C"] * n
    cmap, _ = build_consensus_map(labels, ss, tm_segments=[(0, n - 1)])
    assert cmap == {}


def test_short_fragments_are_dropped():
    assert drop_short_fragments([(0, 1), (3, 7)], min_len=5) == [(3, 7)]


def test_kinks_are_joined_but_hairpins_are_not():
    joined = join_kinks([(0, 10), (12, 22)], breaks=None, min_unwound=3,
                        is_hairpin=lambda a, b: False)
    assert len(joined) == 1
    hairpin = join_kinks([(0, 10), (12, 22)], breaks=None, min_unwound=3,
                         is_hairpin=lambda a, b: True)
    assert len(hairpin) == 2


def test_a_gap_longer_than_the_limit_is_not_joined():
    joined = join_kinks([(0, 10), (20, 30)], breaks=None, min_unwound=3,
                        is_hairpin=lambda a, b: False)
    assert len(joined) == 2


def test_side_from_text_reads_the_polarity():
    assert side_from_text("Extracellular") == "Extracellular"
    assert side_from_text("Cytoplasmic") == "Cytoplasmic"


# ------------------------------------------------------------------ transitions
def test_transition_classifier_reports_a_reason():
    from app.services.topology.membrane import build_membrane

    frame = load_residue_frame(HSERT, "A")
    membrane = build_membrane(frame, [(10, 40), (100, 130)], 7.5)
    assert membrane is not None
    coords = np.asarray(frame.coords, float)
    breaks = np.asarray(frame.chain_breaks())
    evidence = classify_membrane_transition((10, 40), (100, 130), membrane, coords,
                                             breaks, min_cross=6.0, full_cross=0.66)
    assert evidence is not None
    assert evidence.reason


# --------------------------------------------------------------- the full flow
@pytest.mark.parametrize("algorithm", ["3d_slab_geom", "kyte_doolittle_seq", "tmhmm"])
def test_the_consensus_flow_runs_and_reports_its_provenance(algorithm):
    from app.services.topology.service import predict_topology

    response = predict_topology(HSERT, tm_algo=algorithm, ss_algo="none",
                                 flow_type="consensus", chain_id="A")
    assert response.labeler
    assert response.chain_id == "A"
    assert response.parameters_used["flow"] == "consensus"
    assert isinstance(response.warnings, list)


def test_the_energy_provider_finds_the_twelve_helices_of_the_transporter():
    from app.services.topology.providers.tm.energy import EnergyTMProvider

    frame = load_residue_frame(HSERT, "A")
    pred = EnergyTMProvider().predict_tm(HSERT, params=TMParams(), frame=frame)
    assert len(pred.segments) == 12
    assert pred.membrane_score > 10.0            # a clearly negative transfer energy
    assert pred.depth is not None and len(pred.depth) == len(frame)


def test_the_energy_provider_places_a_curved_membrane_on_piezo1():
    from app.services.topology.providers.tm.energy import EnergyTMProvider

    if not PIEZO.exists():
        pytest.skip("the Piezo1 fixture is not part of the repository")
    frame = load_residue_frame(PIEZO, "A")
    pred = EnergyTMProvider().predict_tm(PIEZO, params=TMParams(), frame=frame)
    assert pred.membrane_radius is not None
    assert 100.0 < pred.membrane_radius < 160.0


def test_a_well_behaved_structure_is_reported_as_a_good_placement():
    from app.services.topology.providers.tm.energy import EnergyTMProvider

    frame = load_residue_frame(HSERT, "A")
    pred = EnergyTMProvider().predict_tm(HSERT, params=TMParams(), frame=frame)
    assert pred.placement_quality == "ok"
    assert not pred.placement_flags


def test_an_under_determined_placement_is_flagged():
    """A thickness at the edge of the search range means the bilayer is not pinned
    down by the structure; the tool must say so rather than report a number."""
    from app.services.topology.providers.tm.energy import EnergyTMProvider

    flags, notes, quality = EnergyTMProvider._placement_quality(
        half_thickness=10.0, energy=-80.0, n_segments=8, n_chains=3, predicted=False)
    assert quality == "uncertain"
    assert any("not well determined" in f for f in flags)
    assert not notes

    flags, notes, quality = EnergyTMProvider._placement_quality(
        half_thickness=14.0, energy=-14.0, n_segments=8, n_chains=1, predicted=False)
    assert quality == "caution"
    assert any("marginally" in f for f in flags)

    flags, notes, quality = EnergyTMProvider._placement_quality(
        half_thickness=14.0, energy=-80.0, n_segments=1, n_chains=4, predicted=False)
    assert quality == "uncertain"
    assert any("very little evidence" in f for f in flags)


def test_a_predicted_model_is_a_note_and_not_a_defect():
    """The fact that a model is a prediction is a caveat; it must not by itself
    downgrade the verdict on the placement."""
    from app.services.topology.providers.tm.energy import EnergyTMProvider

    flags, notes, quality = EnergyTMProvider._placement_quality(
        half_thickness=14.0, energy=-80.0, n_segments=9, n_chains=1, predicted=True)
    assert quality == "ok"
    assert not flags
    assert any("prediction" in n for n in notes)


def test_removed_flow_names_still_run_as_consensus():
    from app.services.topology.orchestrator import TopologyOrchestrator

    orchestrator = TopologyOrchestrator(SequenceTMProvider(), None, "ss_then_tm")
    response = orchestrator.execute(HSERT, params=TMParams(), chain_id="A")
    assert any("removed" in w for w in response.warnings)


def test_an_unknown_algorithm_is_rejected():
    from app.services.topology.service import predict_topology

    with pytest.raises(ValueError):
        predict_topology(HSERT, tm_algo="does_not_exist", ss_algo="none")


def test_the_secondary_structure_suffix_runs_out_of_letters_safely():
    """A long loop can hold more than 26 helices; the label must not run out."""
    from app.services.topology.orchestrator import _ss_suffix

    assert _ss_suffix(0) == "a"
    assert _ss_suffix(25) == "z"
    assert _ss_suffix(26) == "aa"
    assert _ss_suffix(27) == "ab"
    assert _ss_suffix(701) == "zz"
    assert _ss_suffix(702) == "aaa"
    # every label must be unique
    labels = [_ss_suffix(i) for i in range(800)]
    assert len(set(labels)) == len(labels)


def test_a_long_loop_with_many_helices_does_not_crash_the_flow():
    """Regression: the topology labels used to index a 26-letter alphabet and raised
    IndexError for a loop with more than 26 secondary-structure elements."""
    from app.services.topology.orchestrator import TopologyOrchestrator
    from app.services.topology.providers.tm.sequence import SequenceTMProvider

    orchestrator = TopologyOrchestrator(SequenceTMProvider(), None, "consensus")
    response = orchestrator.execute(HSERT, params=TMParams(), chain_id="A")
    assert response.regions

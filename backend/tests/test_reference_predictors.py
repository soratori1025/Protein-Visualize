"""
Tests for the reference predictors used in the comparison.

    python -m pytest backend/tests -q

These are re-implementations of published rules, so what has to be tested is that the
rules are implemented as published: the consecutive-window requirement of TopPred, the
dissimilarity filter and the majority vote of the TOPCONS-style consensus, and the
correctness of the hydropathy scale the window rule uses.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
for candidate in (ROOT / "bench", ROOT / "backend", ROOT):
    if str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

from app.core.constants import KYTE_DOOLITTLE  # noqa: E402
from reference_predictors import (  # noqa: E402
    CONSENSUS_MIN_LENGTH, TOPPRED_THRESHOLD, TOPPRED_WINDOW,
    _segments_agree, consensus_segments, hydropathy_profile, toppred_segments,
)


def test_the_hydropathy_scale_is_the_published_kyte_doolittle_table():
    assert KYTE_DOOLITTLE["ILE"] == pytest.approx(4.5)
    assert KYTE_DOOLITTLE["LEU"] == pytest.approx(3.8)
    assert KYTE_DOOLITTLE["VAL"] == pytest.approx(4.2)
    assert KYTE_DOOLITTLE["ARG"] == pytest.approx(-4.5)
    assert KYTE_DOOLITTLE["LYS"] == pytest.approx(-3.9)


def test_the_window_profile_uses_the_three_letter_scale():
    # a run of leucines must give a strongly positive window mean
    profile = hydropathy_profile("L" * 19)
    assert profile[9] == pytest.approx(3.8, abs=0.01)
    # and a run of arginines a strongly negative one
    profile = hydropathy_profile("R" * 19)
    assert profile[9] == pytest.approx(-4.5, abs=0.01)


def test_the_profile_is_zero_outside_the_chain():
    profile = hydropathy_profile("A" * 5)
    assert len(profile) == 5
    assert all(0.0 <= p <= 4.5 for p in profile)


def test_toppred_needs_two_consecutive_favourable_windows():
    # a single hydrophobic stretch that only produces one favourable window
    # must not be called a transmembrane segment
    sequence = "R" * 20 + "L" * 3 + "R" * 20
    assert toppred_segments(sequence) == []


def test_toppred_finds_a_hydrophobic_stretch():
    sequence = "R" * 25 + "L" * 25 + "R" * 25
    segments = toppred_segments(sequence)
    assert len(segments) == 1
    start, end = segments[0]
    assert start <= 25 and end >= 48


def test_toppred_merges_two_plateaus_separated_by_one_dip():
    # two hydrophobic stretches separated by a short hydrophilic stretch must
    # give one segment, not two overlapping ones
    sequence = "R" * 20 + "L" * 22 + "R" * 2 + "L" * 22 + "R" * 20
    segments = toppred_segments(sequence)
    assert len(segments) == 1
    for a, b in segments:
        for c, d in segments[1:]:
            assert b < c or d < a          # no overlap


def test_two_segments_agree_only_when_they_overlap_enough():
    assert _segments_agree((0, 20), (5, 25))
    assert _segments_agree((0, 10), (5, 15))
    assert not _segments_agree((0, 10), (40, 50))
    assert not _segments_agree((0, 10), (11, 20))


def test_the_consensus_drops_a_method_that_agrees_with_nobody():
    providers = {
        "a": [(0, 20), (30, 50)],
        "b": [(0, 20), (30, 50)],
        "outlier": [(500, 520), (600, 620)],
    }
    segments = consensus_segments("A" * 700, providers, 700)
    # the outlier is filtered out, so the two agreeing methods decide
    assert segments
    assert all(s < 100 for s, _ in segments)


def test_the_consensus_takes_the_majority():
    providers = {
        "a": [(0, 20)],
        "b": [(0, 20)],
        "c": [(100, 120)],
    }
    segments = consensus_segments("A" * 200, providers, 200)
    assert segments == [(0, 20)]


def test_the_consensus_drops_segments_below_the_minimum_length():
    providers = {"a": [(0, 5)], "b": [(0, 5)]}
    assert consensus_segments("A" * 200, providers, 200) == []
    providers = {"a": [(0, CONSENSUS_MIN_LENGTH)], "b": [(0, CONSENSUS_MIN_LENGTH)]}
    assert consensus_segments("A" * 200, providers, 200) == [(0, CONSENSUS_MIN_LENGTH)]


def test_the_constants_are_the_published_ones():
    assert TOPPRED_WINDOW == 19
    assert TOPPRED_THRESHOLD == pytest.approx(1.6)

"""
Viterbi / forward / backward for the TMHMM 2.0 model - numpy port of pyTMHMM's
Cython module (pyTMHMM 1.3.6, pyTMHMM/hmm.pyx, MIT licence, see LICENSE).

Why a port instead of ``pip install pyTMHMM``
---------------------------------------------
pyTMHMM ships a Cython extension and, on PyPI, a single macOS-arm64 wheel. Everywhere
else it is compiled at install time, which needs a C compiler (MSVC Build Tools on
Windows) and fails with NumPy >= 2 (``np.int_t`` no longer exists in NumPy's Cython
declarations; the original tmhmm.py also uses ``np.int``, removed in NumPy 1.24).

The port keeps the reference algorithm step for step; only the inner loop over the
previous state ``k`` is vectorised. Element-wise it performs the same double-precision
additions in the same order, and ``argmax`` keeps the FIRST maximum exactly like the
Cython ``if prob > max_state_prob`` - so the Viterbi path is bit-identical
(checked against the compiled extension, see scripts/check_tmhmm_port.py).

Two upstream defects are corrected (``reference_quirks=True`` reproduces both, for
that check):

1. Viterbi end state. The Cython code keeps two DP rows, the last position is
   written to row ``(L-1) % 2``, but the end state is taken as the argmax of row
   ``L % 2`` - the scores of position L-2. The path stays a valid path, it is just not
   always the most probable one near the C terminus.
2. Posterior. The Cython "backward" multiplies by the emission of the CURRENT position
   (``M[i, j] = sum_k M[i+1, k] T[j, k] * e_j(x_i)``) instead of the next one, so
   forward*backward counts e_j(x_i) twice and the rows no longer sum to one (upstream
   renormalises them). Here the textbook scaled backward is used:
   beta_i(j) = sum_k T[j, k] e_k(x_{i+1}) beta_{i+1}(k) / c_{i+1}; then
   alpha_hat * beta_hat is the exact posterior and every row sums to 1.

Additions (upstream raises KeyError on anything outside the 20 letters):
unknown residues ``X`` (and gaps the caller fills with ``X`` to keep the spacing of
an unresolved loop) get emission 1 in every state, i.e. they are marginalised - the
model still sees the right number of positions but no residue identity. Ambiguity
codes use the sum of the letters they stand for (B = D|N, Z = E|Q, J = I|L) and the
rare amino acids their parent (U -> C, O -> K).
"""
from __future__ import annotations

import numpy as np

AMBIGUOUS = {"B": "DN", "Z": "EQ", "J": "IL", "U": "C", "O": "K"}


def emission_table(sequence: str, emissions: np.ndarray, char_map: dict) -> np.ndarray:
    """(L, S) emission probability of every residue in every state."""
    n_states = emissions.shape[0]
    columns = {}
    out = np.empty((len(sequence), n_states))
    for i, ch in enumerate(sequence.upper()):
        col = columns.get(ch)
        if col is None:
            if ch in char_map:
                col = emissions[:, char_map[ch]]
            elif ch in AMBIGUOUS:
                col = emissions[:, [char_map[c] for c in AMBIGUOUS[ch]]].sum(axis=1)
            else:                           # X, '-', '*', anything unknown
                col = np.ones(n_states)
            columns[ch] = col
        out[i] = col
    return out


def viterbi(emit: np.ndarray, initial: np.ndarray, transitions: np.ndarray,
            reference_quirks: bool = False) -> np.ndarray:
    """Most probable state path (state indices), log space."""
    n_obs, n_states = emit.shape
    with np.errstate(divide="ignore"):
        log_init = np.log(initial)
        log_trans = np.log(transitions)
        log_emit = np.log(emit)
    back = np.zeros((n_obs, n_states), dtype=np.intp)
    cols = np.arange(n_states)
    prev = log_init + log_emit[0]
    before_last = np.zeros(n_states)        # upstream reads this row for the end state
    for i in range(1, n_obs):
        cand = prev[:, None] + log_trans    # cand[k, j] = M[i-1, k] + T[k, j]
        best = np.argmax(cand, axis=0)      # first maximum, as the Cython strict '>'
        back[i] = best
        before_last = prev
        prev = cand[best, cols] + log_emit[i]
    if reference_quirks:
        # upstream: argmax(M[L % 2]) = row of position L-2 (zeros when L == 1)
        state = int(np.argmax(before_last if n_obs > 1 else np.zeros(n_states)))
    else:
        state = int(np.argmax(prev))
    path = np.empty(n_obs, dtype=np.intp)
    for i in range(n_obs - 1, -1, -1):
        path[i] = state
        state = back[i, state]
    return path


def forward(emit: np.ndarray, initial: np.ndarray, transitions: np.ndarray):
    """Scaled forward table alpha_hat (rows sum to 1) and the scaling constants c."""
    n_obs, n_states = emit.shape
    alpha = np.zeros((n_obs, n_states))
    const = np.zeros(n_obs)
    row = initial * emit[0]
    const[0] = row.sum()
    alpha[0] = row / const[0]
    for i in range(1, n_obs):
        row = (alpha[i - 1] @ transitions) * emit[i]
        const[i] = row.sum()
        alpha[i] = row / const[i]
    return alpha, const


def backward(emit: np.ndarray, transitions: np.ndarray, const: np.ndarray,
             reference_quirks: bool = False) -> np.ndarray:
    """Scaled backward table; alpha_hat * beta_hat is the posterior."""
    n_obs, n_states = emit.shape
    beta = np.zeros((n_obs, n_states))
    if reference_quirks:                    # upstream recursion (see module docstring)
        beta[n_obs - 1] = 1.0 / const[n_obs - 1]
        for i in range(n_obs - 2, -1, -1):
            beta[i] = (transitions @ beta[i + 1]) * emit[i] / const[i]
        return beta
    beta[n_obs - 1] = 1.0
    for i in range(n_obs - 2, -1, -1):
        beta[i] = (transitions @ (emit[i + 1] * beta[i + 1])) / const[i + 1]
    return beta


def posterior(emit: np.ndarray, initial: np.ndarray, transitions: np.ndarray,
              reference_quirks: bool = False) -> np.ndarray:
    """(L, S) state posterior P(state_i = j | sequence)."""
    alpha, const = forward(emit, initial, transitions)
    beta = backward(emit, transitions, const, reference_quirks)
    post = alpha * beta
    if reference_quirks:
        post = post / post.sum(axis=1, keepdims=True)
    return post

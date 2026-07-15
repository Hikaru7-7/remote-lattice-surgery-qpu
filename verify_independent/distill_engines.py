"""Engines for INDEPENDENT verification of Bell-pair distillation circuits.

Written from scratch 2026-07-14 for adversarial re-derivation of the thesis
"Option B" bit double-selection step (3 Bell pairs -> 1). No project code reused.

CONVENTIONS (fixed throughout; handedness of keep rules is derived from these):
  |Phi+> = (|00>+|11>)/sqrt(2)   (nominal / no-error state)
  Error label P in {I,X,Y,Z} means the pair state (I_Alice (x) P_Bob)|Phi+>:
      I=(a=0,b=0) -> Phi+ = (|00>+|11>)/sqrt2
      X=(a=1,b=0) -> Psi+ = (|01>+|10>)/sqrt2      (bit flip)
      Y=(a=1,b=1) -> Psi- = (|01>-|10>)/sqrt2      (bit+phase flip; phase of -i dropped)
      Z=(a=0,b=1) -> Phi- = (|00>-|11>)/sqrt2      (phase flip)
  a = X-part of the label, b = Z-part of the label (P ~ X^a Z^b up to phase).
  Fidelity F = q_I.

  Bilateral CNOT(c -> t): Alice applies CNOT(A_c -> A_t) AND Bob applies
  CNOT(B_c -> B_t). It maps |Phi+>|Phi+> -> |Phi+>|Phi+> and conjugates labels:
      a_t ^= a_c   (bit flips copy control -> target)
      b_c ^= b_t   (phase flips copy target -> control)

  Bilateral Z measurement of pair i (Alice & Bob each measure Z, compare):
      outcomes AGREE  <=>  a_i = 0   (Phi+/Phi- correlated in Z; Psi+/- anti)
  Bilateral X measurement of pair i:
      outcomes AGREE  <=>  b_i = 0   (Phi+/Psi+ correlated in X; Phi-/Psi- anti)

Two independent engines:
  1. run_label_circuit : exact Bell-diagonal probability-vector simulation by
     exhaustive enumeration over Pauli-label patterns (floats/Fractions/sympy).
  2. run_dm_circuit    : exact density-matrix simulation (2n qubits, numpy).
"""
from itertools import product
import numpy as np

LABELS = [(0, 0), (1, 0), (1, 1), (0, 1)]  # order: I, X, Y, Z
LNAME = {(0, 0): 'I', (1, 0): 'X', (1, 1): 'Y', (0, 1): 'Z'}
NAME2LABEL = {v: k for k, v in LNAME.items()}


def dist_from_q(q4):
    """q4 = (qI, qX, qY, qZ) -> dict label -> prob."""
    return {(0, 0): q4[0], (1, 0): q4[1], (1, 1): q4[2], (0, 1): q4[3]}


# ---------------------------------------------------------------------------
# Circuits. Pair indices: 0 = survivor, 1 = "Z-ancilla" pair, 2 = "X-ancilla".
# ---------------------------------------------------------------------------
OPTION_B = dict(
    n=3,
    gates=[('CNOT', 0, 1), ('CNOT', 2, 1)],   # (name, control_pair, target_pair)
    checks=[(1, 'Z'), (2, 'X')],              # measure pair 1 in Z, pair 2 in X
    name='Option B double selection (CNOT S->Zanc, CNOT Xanc->Zanc)')

OPTION_B_SWAPPED_ORDER = dict(
    n=3,
    gates=[('CNOT', 2, 1), ('CNOT', 0, 1)],
    checks=[(1, 'Z'), (2, 'X')],
    name='Option B with the two CNOTs applied in the opposite order')

SINGLE_SELECTION = dict(
    n=2,
    gates=[('CNOT', 0, 1)],
    checks=[(1, 'Z')],
    name='single selection (BBPSSW/DEJMPS core: CNOT S->anc, measure anc in Z)')

WRONG_SECOND_CNOT = dict(
    n=3,
    gates=[('CNOT', 0, 1), ('CNOT', 1, 2)],   # adversarial: Zanc is CONTROL of 2nd CNOT
    checks=[(1, 'Z'), (2, 'X')],
    name='adversarial variant: 2nd CNOT reversed (Zanc->Xanc), same measurements')


# ---------------------------------------------------------------------------
# Engine 1: Bell-diagonal label engine (exact, works with float/Fraction/sympy)
# ---------------------------------------------------------------------------
def propagate(labels, gates):
    """Conjugate a tuple of per-pair labels (a,b) through bilateral CNOTs."""
    lab = [list(l) for l in labels]
    for g, c, t in gates:
        if g != 'CNOT':
            raise ValueError('unknown gate ' + str(g))
        lab[t][0] ^= lab[c][0]   # X copies control -> target
        lab[c][1] ^= lab[t][1]   # Z copies target -> control
    return tuple(tuple(l) for l in lab)


def syndromes(final_labels, checks):
    """Syndrome bit per check: 0 = outcomes agree, 1 = disagree."""
    out = []
    for pair, basis in checks:
        a, b = final_labels[pair]
        out.append(a if basis == 'Z' else b)
    return tuple(out)


def run_label_circuit(dists, gates, checks, survivor=0, keep=None):
    """Exhaustive exact simulation.

    dists: list (len n) of dicts label->prob (float, Fraction, or sympy expr).
    keep : required syndrome tuple; default all zeros (= all checks 'agree').
    Returns dict with:
      p_succ      : total kept probability
      out         : normalized survivor label distribution (dict), or None if p=0
      out_unnorm  : unnormalized survivor distribution in the keep sector
      records     : input pattern -> (syndrome, kept?, survivor_out_label)
      sectors     : syndrome -> {'p': prob, 'out': unnormalized survivor dist}
    """
    n = len(dists)
    if keep is None:
        keep = tuple(0 for _ in checks)
    records = {}
    sectors = {}
    for pattern in product(LABELS, repeat=n):
        w = 1
        for i, lab in enumerate(pattern):
            w = w * dists[i][lab]
        fin = propagate(pattern, gates)
        syn = syndromes(fin, checks)
        outlab = fin[survivor]
        records[pattern] = (syn, syn == keep, outlab)
        sec = sectors.setdefault(syn, {'p': 0, 'out': {l: 0 for l in LABELS}})
        sec['p'] = sec['p'] + w
        sec['out'][outlab] = sec['out'][outlab] + w
    ksec = sectors.get(keep, {'p': 0, 'out': {l: 0 for l in LABELS}})
    p = ksec['p']
    out_unnorm = dict(ksec['out'])
    out = None
    if p != 0:
        out = {l: v / p for l, v in out_unnorm.items()}
    return {'p_succ': p, 'out': out, 'out_unnorm': out_unnorm,
            'records': records, 'sectors': sectors}


# ---------------------------------------------------------------------------
# Engine 2: exact density-matrix engine (numpy). Qubit order:
#   pair i -> Alice qubit 2i, Bob qubit 2i+1 ; qubit 0 is the MSB of the index.
# ---------------------------------------------------------------------------
_SQ2 = 1.0 / np.sqrt(2.0)
BELL_VECS = [np.array([1, 0, 0, 1]) * _SQ2,    # I : Phi+
             np.array([0, 1, 1, 0]) * _SQ2,    # X : Psi+
             np.array([0, 1, -1, 0]) * _SQ2,   # Y : Psi-
             np.array([1, 0, 0, -1]) * _SQ2]   # Z : Phi-
X2 = np.array([[0., 1.], [1., 0.]])
Z2 = np.array([[1., 0.], [0., -1.]])
I2 = np.eye(2)


def bell_diag_dm(q4):
    rho = np.zeros((4, 4))
    for q, v in zip(q4, BELL_VECS):
        rho += float(q) * np.outer(v, v)
    return rho


def cnot_perm(nq, c, t):
    """CNOT(qubit c -> qubit t) as a permutation matrix on 2^nq dims."""
    dim = 1 << nq
    U = np.zeros((dim, dim))
    for i in range(dim):
        j = i ^ (1 << (nq - 1 - t)) if (i >> (nq - 1 - c)) & 1 else i
        U[j, i] = 1.0
    return U


def op_on(nq, ops):
    m = np.array([[1.0]])
    for q in range(nq):
        m = np.kron(m, ops.get(q, I2))
    return m


def run_dm_circuit(qvecs, gates, checks, survivor=0, keep_signs=None):
    """Exact density-matrix run. keep_signs: +1 ('agree') / -1 ('disagree') per
    check; default all +1. Returns (p_kept, [qI', qX', qY', qZ'] of survivor)."""
    n = len(qvecs)
    nq = 2 * n
    dim = 1 << nq
    rho = np.array([[1.0]])
    for q in qvecs:
        rho = np.kron(rho, bell_diag_dm(q))
    for g, c, t in gates:
        assert g == 'CNOT'
        U = cnot_perm(nq, 2 * c, 2 * t) @ cnot_perm(nq, 2 * c + 1, 2 * t + 1)
        rho = U @ rho @ U.T
    if keep_signs is None:
        keep_signs = tuple(1 for _ in checks)
    Pi = np.eye(dim)
    for (pair, basis), s in zip(checks, keep_signs):
        P = X2 if basis == 'X' else Z2
        PP = op_on(nq, {2 * pair: P, 2 * pair + 1: P})
        Pi = Pi @ (np.eye(dim) + s * PP) / 2.0
    sel = Pi @ rho @ Pi.T
    p = float(np.trace(sel))
    if p < 1e-15:
        return p, None
    tns = sel.reshape([2] * (2 * nq))
    row = list(range(nq))
    col = list(range(nq, 2 * nq))
    keepq = [2 * survivor, 2 * survivor + 1]
    for q in range(nq):
        if q not in keepq:
            col[q] = row[q]           # trace out
    out_idx = [row[keepq[0]], row[keepq[1]], col[keepq[0]], col[keepq[1]]]
    red = np.einsum(tns, row + col, out_idx).reshape(4, 4) / p
    comps = [float(v @ red @ v) for v in BELL_VECS]
    return p, comps

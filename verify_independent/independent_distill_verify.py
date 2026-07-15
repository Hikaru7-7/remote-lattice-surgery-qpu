#!/usr/bin/env python3
"""Independent adversarial re-derivation of the thesis 'Option B' bit
double-selection distillation step (Fujii-Yamamoto style, 3 Bell pairs -> 1).

Circuit under test (pair 0 = survivor, 1 = Z-ancilla, 2 = X-ancilla), applied
symmetrically by Alice and Bob:
    CNOT(survivor -> Z-ancilla); CNOT(X-ancilla -> Z-ancilla);
    measure Z-ancilla in Z, X-ancilla in X (both sides, compare over 2-way CC);
    keep survivor iff BOTH comparisons agree.

Everything below is derived from scratch with two independent engines
(Bell-label exhaustive enumeration; exact density-matrix QM) that are
cross-validated against each other and against the published BBPSSW/DEJMPS
single-selection recurrence on Werner states.

Outputs: results.txt (human-readable dump) and results.json (machine-readable).
"""
import json
import sys
import platform
from fractions import Fraction
from itertools import product

import numpy as np
import sympy as sp

from engines import (LABELS, LNAME, NAME2LABEL, OPTION_B, OPTION_B_SWAPPED_ORDER,
                     SINGLE_SELECTION, WRONG_SECOND_CNOT,
                     run_label_circuit, run_dm_circuit, dist_from_q)

OUT_LINES = []
RES = {'meta': {'date': '2026-07-14',
                'python': sys.version.split()[0],
                'platform': platform.platform(),
                'numpy': np.__version__, 'sympy': sp.__version__,
                'rng_seed': 20260714}}


def emit(*parts):
    line = ' '.join(str(p) for p in parts)
    OUT_LINES.append(line)
    print(line)


def hdr(title):
    emit()
    emit('=' * 96)
    emit(title)
    emit('=' * 96)


def pat_name(pattern):
    return '(' + ','.join(LNAME[l] for l in pattern) + ')'


# ===========================================================================
hdr('SECTION 0: CONVENTIONS AND CIRCUIT AS SIMULATED')
emit("""  Nominal Bell state |Phi+> = (|00>+|11>)/sqrt2. Label P means (I x P)|Phi+>:
    I -> Phi+ ; X -> Psi+ ; Y -> Psi- ; Z -> Phi- ;  (a,b) = (X-part, Z-part).
  Bilateral CNOT(c->t) label action:  a_t ^= a_c ;  b_c ^= b_t.
  Bilateral Z-meas of a pair: outcomes agree <=> a = 0 (no bit-flip content).
  Bilateral X-meas of a pair: outcomes agree <=> b = 0 (no phase-flip content).
  Option B: gates CNOT(0->1), CNOT(2->1); measure pair1 in Z, pair2 in X;
  keep iff syndrome = (agree, agree). Pair 0 survivor, 1 Z-anc, 2 X-anc.""")

# ===========================================================================
hdr('SECTION 1: ENGINE VALIDATION (claim 7 first: trust the engine before the circuit)')

rng = np.random.default_rng(RES['meta']['rng_seed'])
val = {}
for circ in (OPTION_B, OPTION_B_SWAPPED_ORDER, SINGLE_SELECTION, WRONG_SECOND_CNOT):
    md = 0.0
    for _ in range(25):
        qs = [rng.dirichlet([3, 1, 1, 1]) for _ in range(circ['n'])]
        r1 = run_label_circuit([dist_from_q(q) for q in qs], circ['gates'], circ['checks'])
        p2, o2 = run_dm_circuit(qs, circ['gates'], circ['checks'])
        d = abs(r1['p_succ'] - p2)
        for lab, o in zip(LABELS, o2):
            d = max(d, abs(r1['out'][lab] - o))
        md = max(md, d)
    val[circ['name']] = md
    emit(f"  [1a] label engine vs exact density matrix, 25 random NON-iid Bell-diagonal inputs,")
    emit(f"       {circ['name']}")
    emit(f"       max |difference| (p_succ and all 4 output components) = {md:.3e}")

# 1b: exhaustive pure-Pauli-pattern check of the label algebra against exact QM
recsB = run_label_circuit([dist_from_q([1, 0, 0, 0])] * 3,
                          OPTION_B['gates'], OPTION_B['checks'])['records']
maxdev = 0.0
for pattern in product(LABELS, repeat=3):
    syn, kept, outlab = recsB[pattern]
    qs = []
    for lab in pattern:
        v = [0, 0, 0, 0]
        v[LABELS.index(lab)] = 1
        qs.append(v)
    signs = tuple(1 if s == 0 else -1 for s in syn)
    p, comps = run_dm_circuit(qs, OPTION_B['gates'], OPTION_B['checks'], keep_signs=signs)
    dev = abs(p - 1.0)
    want = LABELS.index(outlab)
    for k in range(4):
        dev = max(dev, abs(comps[k] - (1.0 if k == want else 0.0)))
    maxdev = max(maxdev, dev)
emit()
emit(f"  [1b] EXHAUSTIVE check: all 64 pure Pauli input patterns pushed through the exact")
emit(f"       density-matrix engine; predicted (syndrome, survivor label) from the label")
emit(f"       algebra reproduced with max deviation = {maxdev:.3e}")
val['exhaustive_pure_pauli_dm'] = maxdev

# 1c: BBPSSW/DEJMPS single selection on Werner inputs vs published closed form
emit()
emit('  [1c] single selection (bilateral CNOT + Z-parity keep) on Werner inputs vs the')
emit('       published BBPSSW recurrence  F\' = [F^2 + ((1-F)/3)^2] / [F^2 + 2F(1-F)/3 + 5((1-F)/3)^2]:')
bb_err = 0.0
for F in (0.6, 0.75, 0.85, 0.95):
    e = (1 - F) / 3
    d = dist_from_q([F, e, e, e])
    r = run_label_circuit([d, d], SINGLE_SELECTION['gates'], SINGLE_SELECTION['checks'])
    p_cf = F * F + 2 * F * (1 - F) / 3 + 5 * ((1 - F) / 3) ** 2
    Fp_cf = (F * F + ((1 - F) / 3) ** 2) / p_cf
    bb_err = max(bb_err, abs(r['p_succ'] - p_cf), abs(r['out'][(0, 0)] - Fp_cf))
    emit(f"       F={F:.2f}: engine p={r['p_succ']:.12f} closed={p_cf:.12f} | "
         f"engine F'={r['out'][(0, 0)]:.12f} closed={Fp_cf:.12f}")
emit(f"       max deviation = {bb_err:.3e}")
val['bbpssw_werner_closed_form'] = bb_err
RES['validation'] = val

# ===========================================================================
hdr('SECTION 2: EXACT SYMBOLIC RESULTS (iid Bell-diagonal inputs (qI,qX,qY,qZ))')

qI, qX, qY, qZ, t = sp.symbols('qI qX qY qZ t', nonnegative=True)
dsym = {(0, 0): qI, (1, 0): qX, (1, 1): qY, (0, 1): qZ}

rB = run_label_circuit([dsym] * 3, OPTION_B['gates'], OPTION_B['checks'])
pB = sp.expand(rB['p_succ'])
hand = sp.expand((qI + qZ) * (qI**2 + qX**2 + qY**2 + qZ**2)
                 + 2 * (qX + qY) * (qI * qX + qY * qZ))
emit('  [2a] Option B success probability, exact polynomial:')
emit('       p_succ =', pB)
emit('       equals hand-derived closed form (qI+qZ)(qI^2+qX^2+qY^2+qZ^2) + 2(qX+qY)(qI*qX+qY*qZ):',
     sp.simplify(pB - hand) == 0)
emit()
emit('  [2b] Unnormalized kept-output weights N_L (survivor label L; q\'_L = N_L / p_succ):')
NB = {}
for L in LABELS:
    NB[LNAME[L]] = sp.expand(rB['out_unnorm'][L])
    emit(f"       N_{LNAME[L]} =", NB[LNAME[L]])

# order-graded series: qP -> t*qP, qI -> 1 - t(qX+qY+qZ); t^k = k-th order in error
sub = {qI: 1 - t * (qX + qY + qZ), qX: t * qX, qY: t * qY, qZ: t * qZ}


def ser(expr, n=3):
    return sp.expand(sp.series(expr.subs(sub, simultaneous=True), t, 0, n).removeO())


def by_order(expr_t, kmax=2):
    po = sp.Poly(sp.expand(expr_t), t)
    return {k: sp.expand(po.coeff_monomial(t**k)) for k in range(kmax + 1)}


emit()
emit('  [2c] Series to 2nd order in input error components (t grades the order):')
series_B = {}
pB_t = ser(pB)
series_B['p_succ'] = str(pB_t)
emit('       p_succ   =', pB_t)
for L in LABELS:
    s = ser(rB['out_unnorm'][L] / pB)
    series_B['q_' + LNAME[L]] = str(s)
    emit(f"       q'_{LNAME[L]}     =", s)
bitc = ser((rB['out_unnorm'][(1, 0)] + rB['out_unnorm'][(1, 1)]) / pB)
series_B['bit_content'] = str(bitc)
emit("       q'_X+q'_Y =", bitc, '   <-- output bit-flip content (claim 2)')

# swapped CNOT order: must be polynomial-identical
rBs = run_label_circuit([dsym] * 3, OPTION_B_SWAPPED_ORDER['gates'], OPTION_B_SWAPPED_ORDER['checks'])
same = sp.simplify(sp.expand(rBs['p_succ']) - pB) == 0
for L in LABELS:
    same = same and sp.simplify(sp.expand(rBs['out_unnorm'][L]) - rB['out_unnorm'][L]) == 0
emit()
emit('  [2d] Swapping the order of the two CNOTs leaves p_succ and every N_L polynomial')
emit('       IDENTICAL (the two CNOTs share only the target and commute):', same)

# single selection symbolic
rS = run_label_circuit([dsym] * 2, SINGLE_SELECTION['gates'], SINGLE_SELECTION['checks'])
pS = sp.expand(rS['p_succ'])
series_S = {'p_succ': str(ser(pS))}
emit()
emit('  [2e] SINGLE selection (no X-ancilla) for comparison (claim 4):')
emit('       p_succ =', pS, ' -> series', ser(pS))
for L in LABELS:
    s = ser(rS['out_unnorm'][L] / pS)
    series_S['q_' + LNAME[L]] = str(s)
    emit(f"       q'_{LNAME[L]} =", s)

# wrong-direction variant symbolic
rW = run_label_circuit([dsym] * 3, WRONG_SECOND_CNOT['gates'], WRONG_SECOND_CNOT['checks'])
pW = sp.expand(rW['p_succ'])
series_W = {'p_succ': str(ser(pW))}
emit()
emit('  [2f] ADVERSARIAL variant: 2nd CNOT reversed (Z-anc control -> X-anc target),')
emit('       same measurements/keep rule (claim 8: verification direction):')
for L in LABELS:
    s = ser(rW['out_unnorm'][L] / pW)
    series_W['q_' + LNAME[L]] = str(s)
    emit(f"       q'_{LNAME[L]} =", s)
RES['exact'] = {'p_succ_poly': str(pB), 'N_polys': NB and {k: str(v) for k, v in NB.items()}}
RES['series'] = {'optionB': series_B, 'single_selection': series_S,
                 'wrong_second_cnot': series_W}

# ===========================================================================
hdr('SECTION 3: EXHAUSTIVE PAULI-PATTERN AUDIT (claims 1, 4, 8)')

# claim 1: survivor bit-component invariance, all 64 patterns
claim1_ok = all(outlab[0] == pattern[0][0] for pattern, (syn, kept, outlab) in recsB.items())
claim1_b_leak = any(outlab[1] != pattern[0][1] for pattern, (syn, kept, outlab) in recsB.items())
emit(f"  [3a] Claim 1 (no bit leak): survivor OUT bit-component == survivor IN bit-component")
emit(f"       for ALL 64 Pauli patterns: {claim1_ok}")
emit(f"       ... but the survivor PHASE component IS modified by some patterns (b_S <- b_S xor b_Zanc): {claim1_b_leak}")

emit()
emit('  [3b] Weight-1 audit (exactly one non-identity label among the 3 pairs).')
emit('       pattern=(survivor,Zanc,Xanc); syndrome=(Zcheck,Xcheck), 0=agree; kept iff (0,0):')
w1 = []
for pos in range(3):
    for P in ('X', 'Y', 'Z'):
        pattern = tuple(NAME2LABEL[P] if i == pos else (0, 0) for i in range(3))
        syn, kept, outlab = recsB[pattern]
        w1.append({'pattern': pat_name(pattern), 'syndrome': str(syn),
                   'kept': kept, 'survivor_out': LNAME[outlab]})
        emit(f"       {pat_name(pattern):10s} syndrome={syn}  kept={str(kept):5s}  survivor_out={LNAME[outlab]}")
emit('       => the ONLY weight-1 event that is kept is the survivor\'s own Z error (designed')
emit('          passthrough, claim 3). Every other weight-1 error, on ANY pair, is rejected.')
emit('          In particular an X error on the X-ancilla does NOT escape: it propagates into')
emit('          the Z-ancilla via CNOT(Xanc->Zanc) and flips the Z-parity check (harmless to')
emit('          the survivor, costs yield only).')

emit()
emit('  [3c] Weight-2 audit: all kept patterns with exactly two non-identity labels:')
w2 = []
for pattern, (syn, kept, outlab) in recsB.items():
    wt = sum(1 for l in pattern if l != (0, 0))
    if wt == 2 and kept:
        w2.append({'pattern': pat_name(pattern), 'survivor_out': LNAME[outlab],
                   'corrupts': outlab != (0, 0)})
        emit(f"       {pat_name(pattern):10s} kept, survivor_out={LNAME[outlab]}"
             f"{'   <-- corrupts survivor' if outlab != (0,0) else '   (harmless masking)'}")
emit('       => second-order error mechanisms are exactly these coincidences; nothing else.')
RES['audit'] = {'claim1_survivor_bit_invariant': claim1_ok,
                'survivor_phase_backaction_exists': claim1_b_leak,
                'weight1': w1, 'weight2_kept': w2}

# ===========================================================================
hdr('SECTION 4: NUMERIC TABLES (exact rational arithmetic) - claims 5 and 6')


def werner(F):
    e = (1 - F) / 3
    return [F, e, e, e]


def bitdom(F):
    return [F, 1 - F, 0, 0]


def phasedom(F):
    return [F, 0, 0, 1 - F]


MODELS = {'werner': werner, 'bit-dominant': bitdom, 'phase-dominant': phasedom}
FGRID = [Fraction(90, 100), Fraction(95, 100), Fraction(98, 100), Fraction(99, 100)]
tables = {}
for mname, mfn in MODELS.items():
    tables[mname] = {}
    emit()
    emit(f"  input model: {mname}" +
         {'werner': '   q = (F, e, e, e), e=(1-F)/3',
          'bit-dominant': '   q = (F, 1-F, 0, 0)',
          'phase-dominant': '   q = (F, 0, 0, 1-F)'}[mname])
    cols = ['F_raw', 'p_succ', '3/p_succ', "F_out=qI'", "qX'", "qY'", "qZ'", "bit qX'+qY'"]
    emit(f"    {cols[0]:>6s} {cols[1]:>12s} {cols[2]:>10s} {cols[3]:>12s} "
         f"{cols[4]:>12s} {cols[5]:>12s} {cols[6]:>12s} {cols[7]:>16s}")
    for F in FGRID:
        d = dist_from_q(mfn(F))
        r = run_label_circuit([d] * 3, OPTION_B['gates'], OPTION_B['checks'])
        p = r['p_succ']
        o = r['out']
        row = {'p_succ': float(p), 'rate_3_over_p': float(3 / p),
               'qI': float(o[(0, 0)]), 'qX': float(o[(1, 0)]),
               'qY': float(o[(1, 1)]), 'qZ': float(o[(0, 1)]),
               'bit_content': float(o[(1, 0)] + o[(1, 1)])}
        tables[mname][str(float(F))] = row
        emit(f"    {float(F):6.2f} {row['p_succ']:12.8f} {row['rate_3_over_p']:10.4f} "
             f"{row['qI']:12.8f} {row['qX']:12.3e} {row['qY']:12.3e} {row['qZ']:12.3e} "
             f"{row['bit_content']:16.3e}")
RES['tables'] = tables

# what reproduces the thesis x3.4 rate factor?
emit()
emit('  [4b] Solving 3/p_succ = 3.4 for each input model (bisection):')
rate34 = {}
for mname, mfn in MODELS.items():
    def rate(F, mfn=mfn):
        d = dist_from_q([float(v) for v in mfn(F)])
        return 3.0 / run_label_circuit([d] * 3, OPTION_B['gates'], OPTION_B['checks'])['p_succ']
    lo, hi = 0.751, 0.99999
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if rate(mid) > 3.4:
            lo = mid
        else:
            hi = mid
    Fstar = 0.5 * (lo + hi)
    rate34[mname] = {'F_star': Fstar, 'error_rate': 1 - Fstar}
    emit(f"       {mname:15s}: 3/p = 3.4  at  F_raw = {Fstar:.4f}  (raw error = {1-Fstar:.4f})")
RES['rate34'] = rate34

# ===========================================================================
hdr('SECTION 5: CLAIM 4 DEMONSTRATION - the point of the SECOND selection')

zeta = Fraction(5, 100)
dS = dist_from_q([1, 0, 0, 0])
dZ = dist_from_q([1 - zeta, 0, 0, zeta])
dX = dist_from_q([1, 0, 0, 0])
demo = {}
r1 = run_label_circuit([dS, dZ], SINGLE_SELECTION['gates'], SINGLE_SELECTION['checks'])
r2 = run_label_circuit([dS, dZ, dX], OPTION_B['gates'], OPTION_B['checks'])
r3 = run_label_circuit([dS, dZ, dX], WRONG_SECOND_CNOT['gates'], WRONG_SECOND_CNOT['checks'])
emit(f"  Inputs: survivor PERFECT, X-ancilla PERFECT, Z-ancilla has phase error qZ = {float(zeta)}.")
emit(f"  (The only error present is the one that back-propagates onto the survivor control.)")
emit(f"    single selection      : p_succ = {float(r1['p_succ']):.6f}   q'_Z(out) = {float(r1['out'][(0,1)]):.6f}   <- leaks at 1st order")
emit(f"    Option B (double sel.): p_succ = {float(r2['p_succ']):.6f}   q'_Z(out) = {float(r2['out'][(0,1)]):.6f}   <- REJECTED exactly")
emit(f"    wrong 2nd-CNOT variant: p_succ = {float(r3['p_succ']):.6f}   q'_Z(out) = {float(r3['out'][(0,1)]):.6f}   <- leaks: direction matters")
demo['single'] = {'p': float(r1['p_succ']), 'qZ_out': float(r1['out'][(0, 1)])}
demo['optionB'] = {'p': float(r2['p_succ']), 'qZ_out': float(r2['out'][(0, 1)])}
demo['wrong_variant'] = {'p': float(r3['p_succ']), 'qZ_out': float(r3['out'][(0, 1)])}
emit()
emit('  First-order coefficient of q\'_Z (from Section 2 series):')
emit('    single selection : 2*qZ   (survivor qZ + Z-ancilla qZ both pass)')
emit('    Option B         : 1*qZ   (only the survivor\'s own qZ passes; ancilla qZ needs a')
emit('                       coincident X-ancilla phase error -> 2nd order qZ^2+qY^2 terms)')
RES['claim4_demo'] = demo

# ===========================================================================
hdr('SECTION 6: KEEP-RULE HANDEDNESS vs BELL CONVENTION (claim 7)')

emit('  [6a] Syndrome sectors for Werner F=0.90 inputs (which sector holds the fidelity?):')
dW = dist_from_q([Fraction(9, 10)] + [Fraction(1, 30)] * 3)
rsec = run_label_circuit([dW] * 3, OPTION_B['gates'], OPTION_B['checks'])
hand_tbl = {}
for syn in sorted(rsec['sectors']):
    sec = rsec['sectors'][syn]
    p = sec['p']
    fid = float(sec['out'][(0, 0)] / p)
    hand_tbl[str(syn)] = {'prob': float(p), 'fidelity_if_kept': fid}
    emit(f"       syndrome {syn} (Zcheck,Xcheck; 0=agree): prob = {float(p):.6f}, "
         f"survivor fidelity if kept = {fid:.6f}")
emit('       => (agree, agree) is the correct keep sector for the Phi+ convention.')

emit()
emit('  [6b] Uniform Bell-frame inputs (all three pairs identically the given pure Bell state):')
frame_tbl = {}
for nm in ('I', 'X', 'Y', 'Z'):
    v = [0, 0, 0, 0]
    v[LABELS.index(NAME2LABEL[nm])] = 1
    rf = run_label_circuit([dist_from_q(v)] * 3, OPTION_B['gates'], OPTION_B['checks'], keep=None)
    # find the sector holding all the weight
    syn = [s for s, sec in rf['sectors'].items() if sec['p'] == 1][0]
    outl = [l for l, w in rf['sectors'][syn]['out'].items() if w == 1][0]
    bell = {'I': 'Phi+', 'X': 'Psi+', 'Y': 'Psi-', 'Z': 'Phi-'}
    frame_tbl[bell[nm]] = {'syndrome': str(syn), 'survivor': bell[LNAME[outl]]}
    emit(f"       all pairs {bell[nm]:4s}: syndrome = {syn} "
         f"({'agree' if syn[0]==0 else 'DISAGREE'},{'agree' if syn[1]==0 else 'DISAGREE'})"
         f", survivor = {bell[LNAME[outl]]}")
emit('       => keep-parity handedness DEPENDS on the nominal Bell state: the Z-check keep')
emit('          parity flips to DISAGREE for Psi+/Psi- (singlet-type heralds); the X-check')
emit('          keep parity is agree for all uniform frames (ancilla phases cancel pairwise).')
emit('          A global Phi- frame silently self-corrects to Phi+; a global Psi frame')
emit('          distills to Psi+ under the flipped rule.')
RES['handedness'] = {'werner_sectors_F0.90': hand_tbl, 'uniform_frames': frame_tbl}

# ===========================================================================
hdr('SECTION 7: WRITE-OUT')
with open('/home/claude/independent-verify/distill/results.json', 'w') as f:
    json.dump(RES, f, indent=1)
emit('  wrote results.json')
with open('/home/claude/independent-verify/distill/results.txt', 'w') as f:
    f.write('\n'.join(OUT_LINES) + '\n')
emit('  wrote results.txt')

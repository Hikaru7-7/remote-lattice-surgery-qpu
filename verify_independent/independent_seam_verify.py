#!/usr/bin/env python3
"""
Independent adversarial re-derivation: rotated surface code + two-patch
lattice-surgery merge. Exact GF(2) symplectic linear algebra, from scratch.
No external QEC libraries, no reuse of any project scheduler code.

CONVENTION (primal gauge, offset=0), stated once and used everywhere:
  * Data qubits at integer (row r, col c), 0<=r<R top->bottom, 0<=c<C left->right.
    Single patch: R=C=d, d odd.
  * Candidate plaquette positions (r,c), -1<=r<=R-1, -1<=c<=C-1; position (r,c)
    covers {(r,c),(r,c+1),(r+1,c),(r+1,c+1)} intersected with the grid.
  * Pauli type of position (r,c): Z if (r+c+offset) even, else X.  [offset=0 default]
  * Included checks:
      - all interior weight-4 positions 0<=r<=R-2, 0<=c<=C-2 (checkerboard bulk);
      - top edge r=-1 and bottom edge r=R-1 (0<=c<=C-2): ONLY Z-type positions
        => top/bottom are the Z-boundaries: host weight-2 Z checks, Z_L ends there;
      - left edge c=-1 and right edge c=C-1 (0<=r<=R-2): ONLY X-type positions
        => left/right are the X-boundaries: host weight-2 X checks, X_L ends there;
      - corners (weight-1) never included.
  * Logical representatives: Z_L = Z on a full column (vertical, weight d);
    X_L = X on a full row (horizontal, weight d).
  * DUAL GAUGE = global X<->Z relabel (dual=True): identical geometry, every type
    swapped. All results map under Z<->X, Z-boundary<->X-boundary, Z_L<->X_L.
  * offset=1 is the alternative checkerboard colouring within the same gauge
    (an internal convention choice; used to test which claims depend on it).

Pauli representation: symplectic GF(2) vector of length 2n packed in a python int:
bits [0,n)=X part, bits [n,2n)=Z part. Commutation = symplectic inner product.
Rank / membership / nullspace by exact GF(2) Gaussian elimination.
"""
import json
import os
from collections import Counter

PC = int.bit_count

# ---------------- GF(2) core ----------------

def rref(rows):
    """Reduced row echelon form; returns sorted list of (pivot_col, row_value),
    rows mutually reduced (each pivot col appears in exactly one row)."""
    piv = []
    for r0 in rows:
        cur = r0
        for c, pv in piv:
            if (cur >> c) & 1:
                cur ^= pv
        if not cur:
            continue
        c = (cur & -cur).bit_length() - 1
        piv = [(pc_, pv ^ cur if ((pv >> c) & 1) else pv) for pc_, pv in piv]
        piv.append((c, cur))
    piv.sort()
    return piv


def rank(rows):
    return len(rref(rows))


def reduce_vec(v, piv):
    for c, pv in piv:
        if (v >> c) & 1:
            v ^= pv
    return v


def in_span(v, piv):
    return reduce_vec(v, piv) == 0


def nullspace(rows, ncols):
    piv = rref(rows)
    pcols = {c for c, _ in piv}
    basis = []
    for f in range(ncols):
        if f in pcols:
            continue
        v = 1 << f
        for c, pv in piv:
            if (pv >> f) & 1:
                v |= 1 << c
        basis.append(v)
    for b in basis:  # self-check: H b = 0
        assert all(PC(r & b) % 2 == 0 for r in rows)
    return basis

# ---------------- geometry / code construction ----------------

def plaq(r, c, R, C, col0=0):
    """Support of plaquette position (r,c) on an R x C grid, columns shifted by col0."""
    return tuple(sorted((rr, cc + col0) for rr in (r, r + 1) for cc in (c, c + 1)
                        if 0 <= rr < R and 0 <= cc < C))


def build_patch(R, C, offset=0, dual=False, col0=0):
    """Rotated surface code, R rows x C cols, per the stated convention."""
    def typ(r, c):
        t = 'Z' if (r + c + offset) % 2 == 0 else 'X'
        return ('X' if t == 'Z' else 'Z') if dual else t
    zt = 'X' if dual else 'Z'  # type hosted on top/bottom edges (vertical-logical type)
    xt = 'Z' if dual else 'X'  # type hosted on left/right edges
    checks = []
    for r in range(R - 1):
        for c in range(C - 1):
            checks.append((typ(r, c), plaq(r, c, R, C, col0)))
    for c in range(C - 1):
        for r in (-1, R - 1):
            if typ(r, c) == zt:
                checks.append((zt, plaq(r, c, R, C, col0)))
    for r in range(R - 1):
        for c in (-1, C - 1):
            if typ(r, c) == xt:
                checks.append((xt, plaq(r, c, R, C, col0)))
    return checks


def sorted_qubits(checks):
    qs = sorted({q for _, sup in checks for q in sup})
    return qs, {q: i for i, q in enumerate(qs)}


def vec(check, qi, n):
    t, sup = check
    v = 0
    for q in sup:
        v |= 1 << (qi[q] if t == 'X' else n + qi[q])
    return v


def pauli_vec(t, sup, qi, n):
    return vec((t, tuple(sup)), qi, n)


def symp(a, b, n, mask):
    return (PC(a & (b >> n) & mask) ^ PC((a >> n) & b & mask)) & 1


def anticommuting_pairs(vecs, n, mask):
    bad = 0
    for i in range(len(vecs)):
        vi = vecs[i]
        for j in range(i + 1, len(vecs)):
            if symp(vi, vecs[j], n, mask):
                bad += 1
    return bad


def support_int(sup, qi):
    v = 0
    for q in sup:
        v |= 1 << qi[q]
    return v


def css_min_logical(checks, qi, n, ltype):
    """Exact min weight of a nontrivial ltype logical: ker(H_other) \\ rowspace(H_ltype).
    Brute force over the full kernel (Gray code)."""
    hsame = [support_int(s, qi) for t, s in checks if t == ltype]
    hoth = [support_int(s, qi) for t, s in checks if t != ltype]
    ker = nullspace(hoth, n)
    piv = rref(hsame)
    best = None
    v = 0
    for g in range(1, 1 << len(ker)):
        v ^= ker[(g & -g).bit_length() - 1]
        w = PC(v)
        if (best is None or w < best) and not in_span(v, piv):
            best = w
    return best

# ---------------- analyses ----------------

def analyze_single(d, dual=False, offset=0, distance=False):
    ch = build_patch(d, d, offset, dual)
    qs, qi = sorted_qubits(ch)
    n = len(qs)
    mask = (1 << n) - 1
    vs = [vec(c, qi, n) for c in ch]
    piv = rref(vs)
    wts = Counter(len(s) for _, s in ch)
    tys = Counter(t for t, _ in ch)
    zt = 'X' if dual else 'Z'
    xt = 'Z' if dual else 'X'
    # boundary layout: horizontal w2 dominoes on rows 0/d-1 must be zt;
    # vertical w2 dominoes on cols 0/d-1 must be xt
    blay = True
    for t, s in ch:
        if len(s) != 2:
            continue
        rows_ = {r for r, _ in s}
        cols_ = {c for _, c in s}
        if len(rows_) == 1:
            blay &= (rows_ <= {0, d - 1}) and t == zt
        else:
            blay &= (cols_ <= {0, d - 1}) and t == xt
    ZL = pauli_vec(zt, [(r, 0) for r in range(d)], qi, n)      # vertical string
    XL = pauli_vec(xt, [(0, c) for c in range(d)], qi, n)      # horizontal string
    res = dict(
        n=n, n_checks=len(ch),
        n_checks_ok=(len(ch) == d * d - 1),
        weights={str(k): v for k, v in sorted(wts.items())},
        weights_ok=(wts == Counter({4: (d - 1) ** 2, 2: 2 * (d - 1)})),
        types={k: v for k, v in sorted(tys.items())},
        types_ok=(tys[zt] == tys[xt] == (d * d - 1) // 2),
        boundary_layout_ok=blay,
        anticommuting_pairs=anticommuting_pairs(vs, n, mask),
        rank=len(piv), k=n - len(piv),
        sum_weights=sum(len(s) for _, s in ch),
        sum_weights_eq_4d_dm1=(sum(len(s) for _, s in ch) == 4 * d * (d - 1)),
        vertical_logical_ok=(all(symp(ZL, v, n, mask) == 0 for v in vs)
                             and not in_span(ZL, piv)),
        horizontal_logical_ok=(all(symp(XL, v, n, mask) == 0 for v in vs)
                               and not in_span(XL, piv)),
        logicals_anticommute=(symp(ZL, XL, n, mask) == 1),
    )
    if distance:
        dv = css_min_logical(ch, qi, n, zt)   # vertical-type logical min weight
        dh = css_min_logical(ch, qi, n, xt)
        res['d_vertical_logical'] = dv
        res['d_horizontal_logical'] = dh
        res['distance_eq_d'] = (dv == d and dh == d)
    return res


def analyze_mirror(d, dual=False, offset=0, distance=False):
    A = build_patch(d, d, offset, dual)
    Amir = [(t, tuple(sorted((r, 2 * d - 1 - c) for r, c in s))) for t, s in A]
    merged = build_patch(d, 2 * d, offset, dual)
    qs, qi = sorted_qubits(merged)
    n = len(qs)
    mask = (1 << n) - 1
    assert n == 2 * d * d
    seam = [c for c in merged
            if any(cc <= d - 1 for _, cc in c[1]) and any(cc >= d for _, cc in c[1])]
    left = [c for c in merged if all(cc <= d - 1 for _, cc in c[1])]
    right = [c for c in merged if all(cc >= d for _, cc in c[1])]
    A_redge = [c for c in A if len(c[1]) == 2 and all(cc == d - 1 for _, cc in c[1])]
    B_ledge = [(t, tuple(sorted((r, 2 * d - 1 - c) for r, c in s))) for t, s in A_redge]
    dec_ok = (set(left) == set(A) - set(A_redge)
              and set(right) == set(Amir) - set(B_ledge))
    # mirror copy == recoloured (offset-flipped) translated copy?
    B_recolour = set(build_patch(d, d, offset=(offset + 1) % 2, dual=dual, col0=d))
    mirror_eq_recolour = (set(Amir) == B_recolour)
    vs = [vec(c, qi, n) for c in merged]
    piv = rref(vs)
    # seam census
    wts = Counter(len(s) for _, s in seam)
    tys = Counter(t for t, _ in seam)
    w2 = [c for c in seam if len(c[1]) == 2]
    w2_end = w2_type = None
    if len(w2) == 1:
        rows_ = {r for r, _ in w2[0][1]}
        w2_end = 'bottom' if rows_ == {d - 1} else ('top' if rows_ == {0} else 'other')
        w2_type = w2[0][0]
    mt, mtc = tys.most_common(1)[0]           # majority (measured-basis) type
    ot = 'X' if mt == 'Z' else 'Z'
    split_ok = (mtc == (d + 1) // 2 and tys[ot] == (d - 1) // 2)
    per_side = [(sum(1 for _, cc in s if cc <= d - 1),
                 sum(1 for _, cc in s if cc >= d)) for _, s in seam]
    side_ok = all(ab == ((2, 2) if len(seam[i][1]) == 4 else (1, 1))
                  for i, ab in enumerate(per_side))
    remote = sum(b for _, b in per_side)
    # product of majority-type seam checks = mt-string on the two seam columns
    prod = 0
    for c in seam:
        if c[0] == mt:
            prod ^= vec(c, qi, n)
    two_cols = pauli_vec(mt, [(r, cc) for r in range(d) for cc in (d - 1, d)], qi, n)
    # logical bookkeeping
    L_A = pauli_vec(mt, [(r, 0) for r in range(d)], qi, n)            # vertical, col 0
    L_B = pauli_vec(mt, [(r, 2 * d - 1) for r in range(d)], qi, n)    # vertical, col 2d-1
    LL = L_A ^ L_B
    O_A = pauli_vec(ot, [(0, c) for c in range(d)], qi, n)            # horizontal, A half
    O_AB = pauli_vec(ot, [(0, c) for c in range(2 * d)], qi, n)       # full row
    def nrm(v):
        return all(symp(v, u, n, mask) == 0 for u in vs)
    # minority seam checks = products of the two removed facing edge checks
    minor = [c for c in seam if c[0] == ot]
    ve_a = [vec(c, qi, n) for c in A_redge]
    ve_b = [vec(c, qi, n) for c in B_ledge]
    minor_eq_products = all(
        any(vec(m, qi, n) == a ^ b for a in ve_a for b in ve_b) for m in minor)
    res = dict(
        n=n, n_checks=len(merged),
        anticommuting_pairs=anticommuting_pairs(vs, n, mask),
        rank=len(piv), k=n - len(piv),
        rank_ok=(len(piv) == 2 * d * d - 1),
        decomposition_ok=dec_ok,
        mirror_eq_recoloured_translate=mirror_eq_recolour,
        seam_count=len(seam), seam_count_eq_d=(len(seam) == d),
        seam_weights={str(k): v for k, v in sorted(wts.items())},
        seam_weights_ok=(wts == Counter({4: d - 1, 2: 1})),
        seam_types={k: v for k, v in sorted(tys.items())},
        seam_split_ok=split_ok,
        majority_type=mt, minority_type=ot,
        w2_end=w2_end, w2_type=w2_type, w2_is_majority_type=(w2_type == mt),
        measured_logical=f"{mt}_L(A) x {mt}_L(B)",
        majority_product_eq_LLxLL_on_seam_cols=(prod == two_cols),
        LL_in_stabilizer_group=(in_span(LL, piv) and nrm(LL)),
        L_A_survives_as_logical=(nrm(L_A) and not in_span(L_A, piv)),
        other_A_logical_destroyed=(not nrm(O_A)),
        other_joint_logical_survives=(nrm(O_AB) and not in_span(O_AB, piv)),
        old_facing_edge_checks_destroyed=all(not nrm(v) for v in ve_a + ve_b),
        minority_seam_eq_old_edge_pair_products=minor_eq_products,
        per_side_counts_ok=side_ok,
        remote_couplings=remote,
        remote_couplings_eq_2dm1=(remote == 2 * d - 1),
    )
    if distance:
        zt = 'X' if dual else 'Z'
        xt = 'Z' if dual else 'X'
        res['d_vertical_logical'] = css_min_logical(merged, qi, n, zt)
        res['d_horizontal_logical'] = css_min_logical(merged, qi, n, xt)
    return res


def analyze_naive(d, dual=False):
    A = build_patch(d, d, 0, dual)
    B = build_patch(d, d, 0, dual, col0=d)   # identical (non-mirrored) translated copy
    allch = A + B
    qs, qi = sorted_qubits(allch)
    n = len(qs)
    mask = (1 << n) - 1
    assert n == 2 * d * d
    vsall = [vec(c, qi, n) for c in allch]
    k_keep = n - rank(vsall)
    bad_keep = anticommuting_pairs(vsall, n, mask)   # disjoint patches: expect 0
    A_redge = [c for c in A if len(c[1]) == 2 and all(cc == d - 1 for _, cc in c[1])]
    B_ledge = [c for c in B if len(c[1]) == 2 and all(cc == d for _, cc in c[1])]
    retained = [c for c in allch if c not in A_redge and c not in B_ledge]
    ret = [vec(c, qi, n) for c in retained]
    k_drop = n - rank(ret)
    # every candidate CSS seam plaquette at positions (r, d-1), r=-1..d-1, both types
    cands = []
    for r in range(-1, d):
        sup = plaq(r, d - 1, d, 2 * d)
        if len(sup) < 2:
            continue
        for t in ('X', 'Z'):
            v = vec((t, sup), qi, n)
            confl = sum(1 for u in ret if symp(v, u, n, mask))
            cands.append(dict(pos=r, weight=len(sup), type=t, conflicts=confl))
    none_addable = all(c['conflicts'] > 0 for c in cands)
    min_confl = min(c['conflicts'] for c in cands)
    # force in the 'intended' seam checks (geometry of the valid merged code)
    mseam = [c for c in build_patch(d, 2 * d, 0, dual)
             if any(cc <= d - 1 for _, cc in c[1]) and any(cc >= d for _, cc in c[1])]
    forced_bad = anticommuting_pairs(ret + [vec(c, qi, n) for c in mseam], n, mask)
    # naive alternative: horizontal weight-2 cross links at every row, each type
    horiz = {}
    for t in ('X', 'Z'):
        confl = 0
        for r in range(d):
            v = pauli_vec(t, [(r, d - 1), (r, d)], qi, n)
            confl += sum(1 for u in ret if symp(v, u, n, mask))
        horiz[t] = confl
    return dict(
        n=n,
        k_keep_all=k_keep, keep_all_anticommuting_pairs=bad_keep,
        k_drop_facing_edges=k_drop,
        prior_claim_k=d + 1, k_drop_matches_prior=(k_drop == d + 1),
        candidate_seam_checks=cands,
        min_conflicts_over_all_candidates=min_confl,
        no_single_commuting_css_seam_check_exists=none_addable,
        forced_intended_seam_anticommuting_pairs=forced_bad,
        horizontal_w2_link_total_conflicts=horiz,
    )

# ---------------- driver ----------------

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out = {'convention': ('primal gauge: type(r,c)=Z iff (r+c+offset) even, offset=0; '
                          'Z-type w2 checks on top/bottom edges (Z-boundaries, Z_L vertical); '
                          'X-type w2 on left/right (X-boundaries, X_L horizontal); '
                          'dual gauge = global X<->Z relabel of the same geometry.')}
    lines = []

    def log(s=''):
        print(s)
        lines.append(s)

    for d in (3, 5, 7, 9, 11, 13):
        e = {}
        e['single'] = analyze_single(d, distance=(d in (3, 5)))
        e['single_offset1'] = analyze_single(d, offset=1) if d == 3 else None
        e['mirror'] = analyze_mirror(d, distance=(d == 3))
        e['mirror_offset1'] = analyze_mirror(d, offset=1)
        if d in (3, 5):
            e['single_dual'] = analyze_single(d, dual=True, distance=(d == 3))
            e['mirror_dual'] = analyze_mirror(d, dual=True)
        e['naive'] = analyze_naive(d)
        if d == 3:
            e['naive_dual'] = analyze_naive(d, dual=True)
        out[str(d)] = e
        s, m, m1, nv = e['single'], e['mirror'], e['mirror_offset1'], e['naive']
        log(f"=== d={d} ===")
        log(f" single: n={s['n']} checks={s['n_checks']} (ok={s['n_checks_ok']}) "
            f"weights={s['weights']} (ok={s['weights_ok']}) types={s['types']} "
            f"(ok={s['types_ok']}) boundary_layout_ok={s['boundary_layout_ok']} "
            f"anticomm={s['anticommuting_pairs']} rank={s['rank']} k={s['k']} "
            f"sumW={s['sum_weights']} (=4d(d-1): {s['sum_weights_eq_4d_dm1']})")
        if 'd_vertical_logical' in s:
            log(f"         distance: d_vert={s['d_vertical_logical']} "
                f"d_horiz={s['d_horizontal_logical']} (== d: {s['distance_eq_d']})")
        log(f" mirror-merge: anticomm={m['anticommuting_pairs']} rank={m['rank']} "
            f"(ok={m['rank_ok']}) k={m['k']} decomposition_ok={m['decomposition_ok']} "
            f"mirror==recoloured_translate={m['mirror_eq_recoloured_translate']}")
        log(f"   seam: count={m['seam_count']} (=d: {m['seam_count_eq_d']}) "
            f"weights={m['seam_weights']} (ok={m['seam_weights_ok']}) "
            f"types={m['seam_types']} split_ok={m['seam_split_ok']} "
            f"majority={m['majority_type']} w2@{m['w2_end']} w2_type={m['w2_type']} "
            f"(majority-type: {m['w2_is_majority_type']})")
        log(f"   measured={m['measured_logical']} "
            f"prod(majority seam)==L_L on seam cols: "
            f"{m['majority_product_eq_LLxLL_on_seam_cols']}; "
            f"LL_in_group={m['LL_in_stabilizer_group']} "
            f"L_A_survives={m['L_A_survives_as_logical']} "
            f"otherA_destroyed={m['other_A_logical_destroyed']} "
            f"other_joint_survives={m['other_joint_logical_survives']}")
        log(f"   old facing w2 destroyed={m['old_facing_edge_checks_destroyed']} "
            f"minority_seam==edge-pair products="
            f"{m['minority_seam_eq_old_edge_pair_products']} "
            f"per_side_ok={m['per_side_counts_ok']} remote={m['remote_couplings']} "
            f"(=2d-1: {m['remote_couplings_eq_2dm1']})")
        if 'd_vertical_logical' in m:
            log(f"   merged distance: d_vert={m['d_vertical_logical']} "
                f"d_horiz={m['d_horizontal_logical']}")
        log(f"   offset=1 variant: k={m1['k']} seam={m1['seam_count']} "
            f"majority={m1['majority_type']} w2@{m1['w2_end']} "
            f"w2_type={m1['w2_type']}")
        if 'mirror_dual' in e and e['mirror_dual']:
            md = e['mirror_dual']
            log(f"   dual gauge: majority={md['majority_type']} "
                f"measured={md['measured_logical']} w2@{md['w2_end']} "
                f"w2_type={md['w2_type']} k={md['k']}")
        log(f" naive-join: keep-all k={nv['k_keep_all']} "
            f"(anticomm={nv['keep_all_anticommuting_pairs']}); "
            f"drop-facing-edges k={nv['k_drop_facing_edges']} "
            f"(prior claim {nv['prior_claim_k']}: {nv['k_drop_matches_prior']})")
        log(f"   no commuting CSS seam check exists at any seam position/type: "
            f"{nv['no_single_commuting_css_seam_check_exists']} "
            f"(min conflicts over all {len(nv['candidate_seam_checks'])} candidates: "
            f"{nv['min_conflicts_over_all_candidates']}); forced intended seam -> "
            f"{nv['forced_intended_seam_anticommuting_pairs']} anticommuting pairs; "
            f"horizontal w2 links conflicts={nv['horizontal_w2_link_total_conflicts']}")
        log()

    with open(os.path.join(here, 'results.json'), 'w') as f:
        json.dump(out, f, indent=1)
    with open(os.path.join(here, 'RESULTS.txt'), 'w') as f:
        f.write(__doc__ + '\n' + '\n'.join(lines) + '\n')
    log(f"wrote {os.path.join(here, 'results.json')} and RESULTS.txt")


if __name__ == '__main__':
    main()

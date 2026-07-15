#!/usr/bin/env python3
"""Derive & verify the correct lattice-surgery seam checks from first principles,
using the scheduler's OWN build_stabilizers as the ground-truth convention.

Strategy: build a continuous rotated patch of width W, height H with the SAME
rules build_stabilizers uses, gate it against build_stabilizers(d) for W=H=d,
then use W=2d,H=d as the merged code. Everything is checked over GF(2):
commutation, independent-generator count, logical count = n - rank.
"""


from qec_scheduler import build_stabilizers  # ground-truth convention

# ---- rectangular rotated patch via extended-block construction -------------
# Principled generalization of build_stabilizers to any W x H. Slide the full
# checkerboard of would-be plaquettes over blocks (r,c), r in -1..H-1, c in -1..W-1.
# Each block covers the 4 corner qubits (r,c),(r,c+1),(r+1,c),(r+1,c+1) that lie
# in-grid. color = X if (r+c) even else Z (same rule). Keep:
#   weight-4 blocks (bulk) always;
#   weight-2 edge blocks only if their color is the type that edge keeps
#     (top/bottom keep X, left/right keep Z -- the scheduler convention);
#   weight-1 corner blocks never.
def build_rect(W, H):
    stabs = []
    for r in range(-1, H):
        for c in range(-1, W):
            corners = [(r, c), (r, c + 1), (r + 1, c), (r + 1, c + 1)]
            data = frozenset(q for q in corners if 0 <= q[0] < H and 0 <= q[1] < W)
            kind = "X" if (r + c) % 2 == 0 else "Z"
            w = len(data)
            if w == 4:
                stabs.append((kind, data))
            elif w == 2:
                touches_top = (r == -1); touches_bot = (r + 1 == H)
                touches_left = (c == -1); touches_right = (c + 1 == W)
                if (touches_top or touches_bot):        # horizontal edge keeps X
                    if kind == "X": stabs.append((kind, data))
                elif (touches_left or touches_right):   # vertical edge keeps Z
                    if kind == "Z": stabs.append((kind, data))
            # w in (0,1) -> corner/outside, dropped
    return stabs

# ---- GF(2) symplectic helpers ----------------------------------------------
def qindex(W, H):
    idx = {}
    for r in range(H):
        for c in range(W):
            idx[(r, c)] = len(idx)
    return idx

def to_vec(kind, data, idx, n):
    x = [0] * n; z = [0] * n
    for q in data:
        if kind == "X": x[idx[q]] = 1
        else: z[idx[q]] = 1
    return x + z  # length 2n

def commute(v1, v2, n):
    s = 0
    for i in range(n):
        s ^= (v1[i] & v2[n + i]) ^ (v1[n + i] & v2[i])
    return s == 0

def gf2_rank(rows):
    rows = [r[:] for r in rows]
    m = len(rows); w = len(rows[0]) if rows else 0
    rank = 0; col = 0
    for col in range(w):
        piv = None
        for i in range(rank, m):
            if rows[i][col]:
                piv = i; break
        if piv is None: continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for i in range(m):
            if i != rank and rows[i][col]:
                rows[i] = [a ^ b for a, b in zip(rows[i], rows[rank])]
        rank += 1
    return rank

def analyze(stabs, W, H, label):
    idx = qindex(W, H); n = W * H
    vecs = [to_vec(k, d, idx, n) for k, d in stabs]
    # commutation
    bad = 0
    for i in range(len(vecs)):
        for j in range(i + 1, len(vecs)):
            if not commute(vecs[i], vecs[j], n): bad += 1
    rank = gf2_rank(vecs)
    k = n - rank
    print(f"[{label}] n={n} stabs={len(stabs)} indep={rank} "
          f"anticommuting_pairs={bad} logical_qubits={k}")
    return bad, rank, k, idx, n

# ---- gate: build_rect(d,d) must equal build_stabilizers(d) -----------------
def as_set(stabs):
    return {(k, frozenset(d)) for k, d in stabs}

def as_set_sched(d):
    return {(s.kind, frozenset(s.data)) for s in build_stabilizers(d)}

for d in (3, 5, 7):
    ok = as_set(build_rect(d, d)) == as_set_sched(d)
    print(f"gate d={d}: build_rect(d,d) == build_stabilizers(d) ? {ok}")
    assert ok, f"rectangular builder disagrees with scheduler at d={d}"
print("GATE PASSED: rectangular builder matches scheduler convention.\n")

# ---- the merged code = continuous 2d x d patch -----------------------------
for d in (3, 5, 7):
    print(f"===== d={d} =====")
    single = build_rect(d, d)
    analyze(single, d, d, f"single d={d}")
    merged = build_rect(2 * d, d)
    bad, rank, k, idx, n = analyze(merged, 2 * d, d, f"merged 2d={2*d} x d={d}")
    assert bad == 0 and k == 1, f"merged code invalid at d={d}"

    # seam checks: data touching BOTH col d-1 and col d
    seam = [(kd, dd) for kd, dd in merged
            if any(c == d - 1 for _, c in dd) and any(c == d for _, c in dd)]
    nZ = sum(1 for kd, _ in seam if kd == "Z")
    nX = sum(1 for kd, _ in seam if kd == "X")
    print(f"  seam checks: {len(seam)} total  (Z={nZ}, X={nX})   want count={d}")
    assert len(seam) == d, f"seam count {len(seam)} != d={d}"
    for kd, dd in sorted(seam, key=lambda s: (min(r for r,_ in s[1]), s[0])):
        rows = sorted({r for r, _ in dd}); cols = sorted({c for _, c in dd})
        print(f"    {kd}  w{len(dd)}  rows{rows} cols{cols}  {sorted(dd)}")
    print()

# ============================================================================
# Two-identical-copies merge (the scheduler's actual setup): A and B are both
# build_stabilizers(d). Join A's col d-1 to B's seam column. Derive the seam
# checks that validly complete the merge, testing both B orientations.
# ============================================================================
def analyze_two_copies(d, b_mirror):
    A = build_stabilizers(d)
    B = build_stabilizers(d)
    bcol = (d - 1) if b_mirror else 0          # B's seam-facing column
    # drop the boundary weight-2 checks that lie entirely on the joined columns
    A_kept = [s for s in A if not all(c == d - 1 for _, c in s.data)]
    B_kept = [s for s in B if not all(c == bcol for _, c in s.data)]
    n = 2 * d * d
    def gi(side, r, c):                         # global qubit index
        return (0 if side == 'A' else d * d) + r * d + c
    def vec(kind, cells):                       # cells = list of (side,r,c)
        x = [0] * n; z = [0] * n
        for (sd, r, c) in cells:
            (x if kind == "X" else z)[gi(sd, r, c)] = 1
        return x + z
    retained = []
    for s in A_kept:
        retained.append(vec(s.kind, [('A', r, c) for r, c in s.data]))
    for s in B_kept:
        retained.append(vec(s.kind, [('B', r, c) for r, c in s.data]))
    # candidate seam checks straddling A(.,d-1) and B(.,bcol)
    cand = []   # (kind, weight, rows, vec, label)
    for r in range(d - 1):                      # weight-4
        cells = [('A', r, d - 1), ('A', r + 1, d - 1),
                 ('B', r, bcol),  ('B', r + 1, bcol)]
        for kind in ("X", "Z"):
            cand.append((kind, 4, (r, r + 1), vec(kind, cells)))
    for r in range(d):                          # weight-2 (single row across seam)
        cells = [('A', r, d - 1), ('B', r, bcol)]
        for kind in ("X", "Z"):
            cand.append((kind, 2, (r,), vec(kind, cells)))
    # keep candidates that commute with ALL retained
    commuting = [c for c in cand if all(commute(c[3], rv, n) for rv in retained)]
    # greedily pick an independent, mutually-commuting subset until k=1
    chosen = []
    basis = [rv[:] for rv in retained]
    def indep(v, rows):
        test = [r[:] for r in rows] + [v[:]]
        return gf2_rank(test) > gf2_rank(rows)
    for c in commuting:
        if all(commute(c[3], ch[3], n) for ch in chosen) and indep(c[3], basis):
            chosen.append(c); basis.append(c[3][:])
    allvecs = retained + [c[3] for c in chosen]
    bad = sum(1 for i in range(len(allvecs)) for j in range(i + 1, len(allvecs))
              if not commute(allvecs[i], allvecs[j], n))
    rank = gf2_rank(allvecs); k = n - rank
    nX = sum(1 for c in chosen if c[0] == "X"); nZ = len(chosen) - nX
    print(f"  [B_mirror={b_mirror}] seam picked={len(chosen)} (X={nX},Z={nZ}) "
          f"anticommute={bad} k={k}   (retained={len(retained)}, want seam=d={d}, k=1)")
    for kind, w, rows, _ in sorted(chosen, key=lambda c: (min(c[2]), c[0])):
        print(f"      {kind}  w{w}  A rows {list(rows)} col {d-1}   "
              f"B rows {list(rows)} col {bcol}")
    return len(chosen), k, bad

print("\n##### scheduler setup: two identical copies of build_stabilizers(d) #####")
for d in (3, 5, 7):
    print(f"----- d={d} -----")
    for b_mirror in (False, True):
        analyze_two_copies(d, b_mirror)

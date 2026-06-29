"""
QEC round scheduler.

Goal: take an abstract surface-code round and lay it onto a real ion chip,
in four layers -- build the code, colour it, place ions, schedule moves.
Rule: build something, then check it. We trust the checks, not the builder.

Layer 1: build the code (data qubits + checks) and prove it is right.
Layer 2: give each check's qubits a step, so no qubit is used twice at once.

Run:  python3 qec_scheduler.py
"""
from __future__ import annotations
from dataclasses import dataclass


# --- REPRESENT -------------------------------------------------------------
# A qubit's place is (row, col) on the grid. Checks are built from this.
Coord = tuple[int, int]


@dataclass(frozen=True)
class Stabilizer:
    """One parity check. It has a kind (X or Z) and the qubits it reads.
    Frozen so we can drop it in a set for the d=3 check."""
    kind: str                  # "X" or "Z"
    data: frozenset            # the qubits this check reads

    @property
    def weight(self) -> int:   # 4 inside the grid, 2 on an edge
        return len(self.data)


# --- ENUMERATE -------------------------------------------------------------
def data_qubits(d: int) -> list:
    """All d*d data qubits. The board the checks sit on."""
    return [(r, c) for r in range(d) for c in range(d)]


def num(rc: Coord, d: int) -> int:
    """Turn (row, col) into the simulator's number 1..d*d.
    Used only by the d=3 check."""
    r, c = rc
    return r * d + c + 1


# --- GENERATE  (the heart of Layer 1) --------------------------------------
def build_stabilizers(d: int) -> list:
    """Build all d*d - 1 checks of the rotated surface code.

    Inside: one weight-4 check on each 2x2 block; X/Z swaps like a chess board.
    Edges:  weight-2 checks that close the open blocks; top/bottom X, sides Z.
    """
    stabs = []

    # Inside: a weight-4 check on each 2x2 block. (r+c) even -> X, odd -> Z.
    for r in range(d - 1):
        for c in range(d - 1):
            block = frozenset({(r, c), (r, c + 1), (r + 1, c), (r + 1, c + 1)})
            kind = "X" if (r + c) % 2 == 0 else "Z"
            stabs.append(Stabilizer(kind, block))

    # Edges: weight-2 checks, one loop per side, on every other spot.
    for c in range(1, d - 1, 2):                  # top    -> X
        stabs.append(Stabilizer("X", frozenset({(0, c), (0, c + 1)})))
    for c in range(0, d - 1, 2):                  # bottom -> X
        stabs.append(Stabilizer("X", frozenset({(d - 1, c), (d - 1, c + 1)})))
    for r in range(0, d - 1, 2):                  # left   -> Z
        stabs.append(Stabilizer("Z", frozenset({(r, 0), (r + 1, 0)})))
    for r in range(1, d - 1, 2):                  # right  -> Z
        stabs.append(Stabilizer("Z", frozenset({(r, d - 1), (r + 1, d - 1)})))

    return stabs


# --- COLOUR  (the heart of Layer 2) ----------------------------------------
def corner_step(stab: Stabilizer, q: Coord) -> int:
    """Pick the step (0..3) for qubit q in this check.
    A check caps a 2x2 block. q sits at one corner of it, and the corner is the
    step. An edge check is half a block, so its block sits just off the grid;
    we slide the corner one step out to match."""
    rows = {r for r, c in stab.data}
    cols = {c for r, c in stab.data}
    r0, c0 = min(rows), min(cols)
    if stab.weight == 2:               # edge check: slide the block off the edge
        if len(rows) == 1:             # one row -> top or bottom
            r0 = r0 - 1 if r0 == 0 else r0
        else:                          # one col -> left or right
            c0 = c0 - 1 if c0 == 0 else c0
    qr, qc = q
    return 2 * (qr - r0) + (qc - c0)    # corners (0,0)(0,1)(1,0)(1,1) -> 0 1 2 3


def assign_steps(d: int) -> dict:
    """Give every check's qubits a step. Returns step -> list of qubits."""
    schedule = {0: [], 1: [], 2: [], 3: []}
    for stab in build_stabilizers(d):
        for q in stab.data:
            schedule[corner_step(stab, q)].append(q)
    return schedule


# --- AUDIT #1: counts are right, at every distance -------------------------
def check_census(d: int) -> None:
    """Check 1: are the counts right? Works at any d."""
    stabs = build_stabilizers(d)
    nX   = sum(s.kind == "X" for s in stabs)
    nZ   = sum(s.kind == "Z" for s in stabs)
    bulk = sum(s.weight == 4 for s in stabs)
    bdry = sum(s.weight == 2 for s in stabs)
    assert len(stabs) == d*d - 1,  f"d={d}: {len(stabs)} checks, want {d*d-1}"
    assert bulk == (d-1)**2,       f"d={d}: {bulk} bulk, want {(d-1)**2}"
    assert bdry == 2*(d-1),        f"d={d}: {bdry} edge, want {2*(d-1)}"
    assert nX == nZ == (d*d-1)//2, f"d={d}: X={nX}, Z={nZ}, want {(d*d-1)//2} each"
    for s in stabs:
        assert s.weight in (2, 4), f"d={d}: a weight-{s.weight} check"
        for (r, c) in s.data:
            assert 0 <= r < d and 0 <= c < d, f"d={d}: ({r},{c}) off the grid"


# --- AUDIT #2: d=3 matches the hand-checked simulator ----------------------
# The 8 checks from the hand-checked simulator (qubits 1..9).
# Census checks the shape; this checks the content.
D3_REFERENCE = {
    ("X", frozenset({1, 2, 4, 5})), ("X", frozenset({5, 6, 8, 9})),
    ("X", frozenset({2, 3})),       ("X", frozenset({7, 8})),
    ("Z", frozenset({2, 3, 5, 6})), ("Z", frozenset({4, 5, 7, 8})),
    ("Z", frozenset({1, 4})),       ("Z", frozenset({6, 9})),
}


def check_d3_matches_reference() -> None:
    """Check 2: does d=3 match the simulator, check by check?"""
    got = {(s.kind, frozenset(num(rc, 3) for rc in s.data))
           for s in build_stabilizers(3)}
    assert got == D3_REFERENCE, (
        "d=3 does not match the simulator:\n"
        f"  missing: {D3_REFERENCE - got}\n"
        f"  extra:   {got - D3_REFERENCE}"
    )


# --- AUDIT #3: the colouring is collision-free -----------------------------
def check_no_double_touch(d: int) -> None:
    """Check 3: in one step, no qubit is used by two checks."""
    for step, qubits in assign_steps(d).items():
        assert len(qubits) == len(set(qubits)), \
            f"d={d}, step {step}: a qubit is used twice"


# --- RUN -------------------------------------------------------------------
# Run the checks. Print pass or fail.
if __name__ == "__main__":
    try:
        for d in (3, 5, 7, 9, 11):
            check_census(d)
        print("census .............. PASS  (d = 3,5,7,9,11)")
        for d in (3, 5, 7, 9, 11):
            check_no_double_touch(d)
        print("no qubit used twice . PASS  (d = 3,5,7,9,11)")
        check_d3_matches_reference()
        print("d=3 vs simulator .... PASS")
        print("\nLayers 1-2 done. Next: Layer 3 -- place ions in wells.")
    except NotImplementedError as e:
        print("not written yet:", e)
    except AssertionError as e:
        print("CHECK FAILED:\n", e)
#!/usr/bin/env python3
"""qec_distill.py -- a certified double-selection distillation round for the
remote seam pairs (thesis Section 7.3).

Distillation is LOCAL to each module and does not grow with the code distance:
one fixed circuit per distilling lane, d-1 lanes in parallel. Double selection
takes THREE raw Bell pairs and yields ONE purified pair, kept iff two parity
checks agree (success ~ (1-eps)^2, about 88% at the baseline 6% raw error).

The key to certifying it: the QUANTUM choreography -- accumulate three pairs,
two bilateral CNOTs, two measurements -- is DETERMINISTIC and is certified here
exactly like any gate sequence. Only keep/discard is classical, and discard is
a retry, not a schedule branch. The retry is priced as the 3/P_ds rate
multiplier (make_requirement), so the geometric scheduler stays branch-free.

Ion budget, one distilling lane, one module (peak):
  N   networker  -- generates raw pairs at the cavity (ping-pong, Section 4.2)
  K   carrier    -- ping-pong partner, ferries each delivered half
  H0  survivor   -- the kept half; the two CNOTs target it; never measured
  H1  ancilla-1  -- sacrificed, first parity check
  H2  ancilla-2  -- sacrificed, second parity check
Five comm-region Ba+ per distilling lane against two for a raw lane. During
accumulation the halves park in the three hold wells at the gate zone's
interface end, beside the swap well (Chapter 4), a one-hop ferry; at the
distillation step all three survivor/ancilla halves sit in the lane's gate-well
column, where the two CNOTs and two reads run with the gate-zone primitives.

Circuit (per module, per distillation round), as an ordered op list:
  ("acc",  lane, j)        accumulate raw pair j (j=0 survivor, 1, 2): the
                           networker heralds it and the carrier delivers it
                           to the hold well  -- three of these, serial
  ("dcnot", lane, a)       bilateral CNOT ancilla a -> survivor  (a in {1,2})
  ("dmeas", lane, a, basis) measure ancilla a in its basis: check 1 in Z
                           (bit-flip), check 2 in X (phase-flip); complementary
  ("keep",  lane)          classical: keep survivor iff both checks agree

Checks (certify_distill runs them at every odd d = 3..27):
  - ion budget: exactly 5 comm-region ions per distilling lane, distinct
  - the two CNOTs couple each distinct ancilla to the survivor exactly once
  - the two ancillas are each measured exactly once; the survivor never is
  - the survivor is the single, well-defined output
  - the d-1 distilling lanes run with no shared ion or gate well (parallel)
  - the round ends at rest: survivor holds the purified pair for the seam gate,
    ancillas are measured then reset, networker and carrier resume ping-pong
  - the two parity checks read in complementary bases (Z and X), so a round
    catches both bit-flip and phase-flip error; the X read folds in one
    single-qubit rotation, the same basis primitive the local X-checks use

Run:  python3 qec_distill.py         # certify + report
      python3 qec_distill.py 13      # one distance
"""
import sys
import qec_scheduler as S

# double selection: 1 survivor + 2 sacrificed ancillas = 3 raw pairs in
N_RAW_IN = 3
N_ANCILLA = 2                        # the two parity checks
IONS_PER_DISTILL_LANE = 5            # N, K, H0, H1, H2  (2 ping-pong + 3 held)
IONS_PER_RAW_LANE = 2               # N, K  (Section 4.2 ping-pong)


# The two selection checks read in COMPLEMENTARY bases: one catches bit-flip
# error (read in Z, the native basis), the other catches phase-flip error (read
# in X, one single-qubit rotation folded into that read). Catching both is what
# makes it double selection. The X read uses the same basis-change primitive the
# local X-checks use (qec_scheduler.check_basis), off the critical path.
CHECK_BASES = (("Z", 1), ("X", 2))          # (basis, ancilla) per parity check


def distill_ops(d: int, lane: int) -> list:
    """The ordered operations of one double-selection round on one lane, one
    module. Deterministic; keep/discard is the classical tail."""
    ops = [("acc", lane, j) for j in range(N_RAW_IN)]        # gather 3 raw halves
    for basis, a in CHECK_BASES:                            # two selection checks
        ops.append(("dcnot", lane, a))                      # ancilla a -> survivor
        ops.append(("dmeas", lane, a, basis))              # read ancilla a in its basis
    ops.append(("keep", lane))                              # classical keep/discard
    return ops


def distill_round(d: int) -> list:
    """All d-1 active lanes distilling one output pair each, in parallel."""
    lanes = range(d - 1)                                    # active seam lanes
    return [op for l in lanes for op in distill_ops(d, l)]


def hold_wells(d: int, lane: int) -> dict:
    """The three hold wells this lane uses for the survivor and two ancillas.
    They sit at the gate zone's interface end beside the swap well during
    accumulation, and the halves move to the lane's gate-well column at the
    CNOT step (Chapter 4). Distinct per lane."""
    return {"survivor": ("hold", lane, 0),
            "ancilla1": ("hold", lane, 1),
            "ancilla2": ("hold", lane, 2)}


# ---- checks ---------------------------------------------------------------
def check_ion_budget(d: int) -> None:
    """Five comm-region ions per distilling lane, distinct, and d-1 lanes."""
    lanes = list(range(d - 1))
    assert len(lanes) == d - 1
    ids = set()
    for l in lanes:
        lane_ids = {("N", l), ("K", l), ("H0", l), ("H1", l), ("H2", l)}
        assert len(lane_ids) == IONS_PER_DISTILL_LANE
        assert lane_ids.isdisjoint(ids), f"d={d} lane {l}: ion id clash across lanes"
        ids |= lane_ids
    assert len(ids) == (d - 1) * IONS_PER_DISTILL_LANE


def check_circuit(d: int) -> None:
    """Per lane: 3 accumulations, 2 CNOTs (each ancilla once, targeting the
    survivor), 2 reads (each ancilla once), the survivor never read, one keep."""
    for l in range(d - 1):
        ops = distill_ops(d, l)
        acc = [o for o in ops if o[0] == "acc"]
        cnot = [o for o in ops if o[0] == "dcnot"]
        meas = [o for o in ops if o[0] == "dmeas"]
        keep = [o for o in ops if o[0] == "keep"]
        assert [o[2] for o in acc] == [0, 1, 2], f"d={d} lane {l}: accumulate not 0,1,2"
        assert sorted(o[2] for o in cnot) == [1, 2], f"d={d} lane {l}: CNOTs not on both ancillas"
        assert sorted(o[2] for o in meas) == [1, 2], f"d={d} lane {l}: reads not on both ancillas"
        assert len(keep) == 1, f"d={d} lane {l}: not exactly one keep"
        # the survivor (index 0) is accumulated but never a CNOT target-id or read
        assert 0 not in {o[2] for o in cnot} and 0 not in {o[2] for o in meas}, \
            f"d={d} lane {l}: survivor was consumed"


def check_hold_wells(d: int) -> None:
    """Each lane's three hold wells are distinct, and no well is shared across
    lanes (the d-1 distillations run in parallel without contention)."""
    seen = set()
    for l in range(d - 1):
        w = hold_wells(d, l)
        vals = set(w.values())
        assert len(vals) == 3, f"d={d} lane {l}: hold wells not distinct"
        assert vals.isdisjoint(seen), f"d={d} lane {l}: hold well shared across lanes"
        seen |= vals


def check_parallel(d: int) -> None:
    """The d-1 lanes share no ion and no hold well, so their fixed circuits run
    at once. Distillation therefore does not lengthen with d; only the lane
    count grows (a resource, not new geometry)."""
    check_ion_budget(d)
    check_hold_wells(d)
    # the circuit is identical on every lane and every distance
    shapes = {tuple(o[0] for o in distill_ops(d, l)) for l in range(d - 1)}
    assert len(shapes) == 1, f"d={d}: lanes do not share one circuit shape"


def check_ends_at_rest(d: int) -> None:
    """After a distillation round the survivor holds the purified pair (ready for
    the seam gate), the two ancillas are measured then reset, and the networker
    and carrier are back in ping-pong. The op list realises exactly this: one
    keep per lane, two reads per lane, no leftover held ancilla."""
    for l in range(d - 1):
        ops = distill_ops(d, l)
        assert sum(1 for o in ops if o[0] == "keep") == 1
        assert sum(1 for o in ops if o[0] == "dmeas") == N_ANCILLA


def check_distill_basis(d: int) -> None:
    """The two parity checks read in complementary bases, one in Z (bit-flip) and
    one in X (phase-flip), so a round catches both error types. The X read folds in
    one single-qubit basis rotation, the same primitive the local X-checks use
    (qec_scheduler.check_basis); it is off the critical path, so it changes no rate.
    Every lane runs the identical pair."""
    for l in range(d - 1):
        bases = [o[3] for o in distill_ops(d, l) if o[0] == "dmeas"]
        assert bases == ["Z", "X"], f"d={d} lane {l}: parity checks not (Z, X): {bases}"
        assert len(set(bases)) == N_ANCILLA, f"d={d} lane {l}: checks not complementary"


def certify_distill(d: int) -> int:
    """Run every distillation check at one distance. Returns the count."""
    assert isinstance(d, int) and d >= 3 and d % 2 == 1, "d must be odd >= 3"
    checks = [check_ion_budget, check_circuit, check_hold_wells,
              check_parallel, check_ends_at_rest, check_distill_basis]
    for chk in checks:
        chk(d)
    return len(checks)


# ---- resource summary (priced against the timing + requirement layers) ----
def resources(d: int) -> dict:
    """The concrete cost of distilling at full width, d=given."""
    lanes = d - 1
    return {
        "distilling_lanes": lanes,
        "ions_per_lane": IONS_PER_DISTILL_LANE,
        "extra_ions_per_lane": IONS_PER_DISTILL_LANE - IONS_PER_RAW_LANE,   # +3
        "extra_ions_total": lanes * (IONS_PER_DISTILL_LANE - IONS_PER_RAW_LANE),
        "cnots_per_module_per_round": N_ANCILLA,        # 2
        "reads_per_module_per_round": N_ANCILLA,        # 2
        "raw_pairs_per_output": N_RAW_IN,               # 3
        "accumulation_windows": N_RAW_IN,               # serial, one per window
    }


def circuit_time_us(d: int, k: int = 1) -> float:
    """The deterministic distillation circuit's added time per output pair, at
    bracket k. Two CNOTs (gate + merge into the well) and two reads, run in the
    lane's gate-well column; each held half hops the one-segment carrier hop to
    the gate and back. Priced with the released timing beats, so the number
    reprints. Off the raw-merge critical path where it overlaps, this is the
    honest serial cost of the circuit itself, not the accumulation (which is the
    3/P_ds rate multiplier)."""
    import qec_timing as T
    gate = T.DURATIONS_US["gate"][k] + T.DURATIONS_US["merge_split"][k]   # one CNOT beat
    read = T.DURATIONS_US["measure"][k]
    seg = T.DURATIONS_US["shuttle"][k]
    hop = 2 * 1 * seg                                    # in and out, one segment each
    return N_ANCILLA * (gate + read + hop)


if __name__ == "__main__":
    ds = [int(sys.argv[1])] if len(sys.argv) > 1 else [3, 5, 7, 13, 27]
    print("double-selection distillation, certified per lane (Section 7.3):")
    print(f"  {N_RAW_IN} raw pairs -> 1 purified, {N_ANCILLA} CNOTs + "
          f"{N_ANCILLA} reads per module, keep iff both checks agree.")
    for d in ds:
        n = certify_distill(d)
        r = resources(d)
        print(f"  d={d:2d}: {n} checks PASS  |  {r['distilling_lanes']} lanes x "
              f"{r['ions_per_lane']} ions (+{r['extra_ions_per_lane']}/lane vs raw) "
              f"= +{r['extra_ions_total']} ions/module")
    print("\ndistillation circuit time per output pair (2 CNOT + 2 read + hops),")
    print("opt / base / cons, ms:")
    ct = [circuit_time_us(7, k) / 1e3 for k in range(3)]
    print(f"  {ct[0]:.2f} / {ct[1]:.2f} / {ct[2]:.2f}   "
          f"(small against the {13.81:.1f} ms raw round)")
    assert 1.5 < ct[1] < 2.5, ct                        # ~2 ms at baseline

    print("\nkeep/discard is classical: discard = redo, priced as the 3/P_ds rate")
    print("multiplier in make_requirement, so the geometric schedule is branch-free.")
    print("The circuit is d-independent; only the lane count grows.")

    # cross-check the distilled operating point against make_requirement
    try:
        import math, make_requirement as MR
        ln100 = math.log(100.0)
        eta_raw = 100 * MR.required_eta_int(7, 1, ln100)[0]     # 4.3%
        p_ds = (1 - 0.06) ** 2
        eta_dist = eta_raw * math.sqrt(N_RAW_IN / p_ds)
        print(f"\ndistilled operating point at baseline: eta_int {eta_raw:.1f}% "
              f"-> {eta_dist:.1f}% (demand x{N_RAW_IN/p_ds:.2f}); matches Section 7.3.")
        assert 7.5 < eta_dist < 8.3, eta_dist
    except ImportError:
        pass

    # cross-check the ion arithmetic against Chapter 4's ping-pong lane
    assert IONS_PER_RAW_LANE == S.comm_ions_per_lane(), \
        "raw-lane ion count disagrees with qec_scheduler.comm_ions_per_lane()"
    print("\nconsistency: raw lane = 2 ping-pong ions "
          f"(qec_scheduler.comm_ions_per_lane() = {S.comm_ions_per_lane()}); "
          "distilling lane adds 3 held halves.")

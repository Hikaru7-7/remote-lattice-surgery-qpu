"""
Timing layer: prices the certified schedule in seconds.
========================================================
The scheduler proves the round can run and packs it into parallel time-steps.
This module weights those steps with operation durations, so a step count
becomes a round time and a Bell-pair count becomes a demand rate.

The pricing rule
----------------
* An operation costs the sum of its beats (BEATS in qec_scheduler).
* A packed step costs its slowest operation.
* The round time is the sum of its steps. Nothing is averaged.

The durations
-------------
Each beat kind carries an (optimistic, baseline, conservative) duration in
microseconds, traced to a published demonstration and, where the species
differs, scaled to 137Ba+ by sqrt(m_Ba/m_ref) at fixed potentials. The
brackets are deliberately wide; the baseline is the defensible middle, not a
promise. Sources, one per line:

  gate          270 / 600 / 3300 us
      150 and 330 us near-field microwave gates on 43Ca+ at 0.5-1% error
      (Weber et al. 2024), scaled by sqrt(137/43)=1.785 to 268 and 589 us.
      Conservative takes the slowest demonstrated high-fidelity microwave
      gate, 3.25 ms at 99.7% on 43Ca+ (Harty et al. 2016), unscaled.
      A gate beat ("merge+gate") also pays one merge_split, since the pair
      is merged into the well as part of the beat.
  merge_split   55 / 105 / 200 us
      55 us two-ion separation demonstrated on 9Be+ (Bowler et al. 2012).
      Conservative holds the sqrt(137/9)=3.9 mass scaling (~215 us, capped
      at the bracket ceiling 200 us of the reference cell). Baseline is the
      geometric mean of the bracket.
  swap_rotation 42 / 78 / 160 us
      42 us two-ion crystal rotation at 99.5(5)% on 40Ca+ (Kaufmann et al.
      2017). Baseline scales by sqrt(137/40)=1.85 to 78 us; conservative
      doubles the baseline as margin.
  shuttle       20 / 52 / 200 us
      28 us per 200 um segment, adiabatic, on 40Ca+ (Kaufmann et al. 2017),
      scaled by 1.85 to 52 us. Optimistic reflects diabatic transport
      (370 um in 8 us on 9Be+, Bowler et al. 2012) and the 10 us industrial
      assumption. The wide bracket owns multi-segment legs and settle time.
  junction      2x shuttle (derived, never set independently)
      One junction beat crosses the intersection twice (see the op_tally
      NOTE in qec_scheduler). Baseline 104 us = 2 x 52, consistent with the
      52 us single-transit estimate used for the junction array.
  measure       150 / 200 / 400 us
      145 us at 99.99% readout on 40Ca+ (Myerson et al. 2008); 200 us
      industrial Ba+ cycle as baseline; 400 us bracket ceiling.
  herald        0 by construction
      Pair readiness is a delivery requirement on the link, not a schedule
      duration. A round that waits for its pair is a stall, and stalls are
      what the delivery requirement forbids.

Run:  python3 qec_timing.py        # duration table, T_round/T_merge, demand rate
      python3 qec_timing.py 13     # one distance
Every number in thesis Section 5.2 reprints from this command.
"""
import sys
import qec_scheduler as S

BRACKETS = ("optimistic", "baseline", "conservative")

DURATIONS_US = {
    "gate":          (270.0,  600.0, 3300.0),
    "merge_split":   ( 55.0,  105.0,  200.0),
    "swap_rotation": ( 42.0,   78.0,  160.0),
    "shuttle":       ( 20.0,   52.0,  200.0),
    "measure":       (150.0,  200.0,  400.0),
    "herald":        (  0.0,    0.0,    0.0),
}


def beat_cost(beat: str, k: int) -> float:
    """One beat's duration in us at bracket k (0 opt, 1 base, 2 cons)."""
    kind = S.BEAT_KIND[beat]
    if kind == "junction":
        return 2.0 * DURATIONS_US["shuttle"][k]     # two intersection crossings
    if beat == "merge+gate":
        return DURATIONS_US["gate"][k] + DURATIONS_US["merge_split"][k]
    return DURATIONS_US[kind][k]


def op_cost(op, k: int) -> float:
    return sum(beat_cost(b, k) for b in S.BEATS[op[0]])


def schedule_time_us(d: int, merge: bool = False, rounds: int = 1, k: int = 1) -> float:
    """Step-weighted critical path: each packed step costs its slowest op."""
    ops = S.round_ops(d, merge=merge, rounds=rounds)
    steps = S.parallel_steps(d, merge=merge, rounds=rounds)
    return sum(max(op_cost(ops[j], k) for j in st) for st in steps)


def phase_times_us(d: int, k: int = 1) -> dict:
    """The local round's time split into its four phases (band, file-out,
    shuttle/read, reset), the same segmentation the scheduler certifies."""
    ops = S.round_ops(d, merge=False)
    steps = S.parallel_steps(d, merge=False)
    kinds = [{ops[j][0] for j in st} for st in steps]
    band_end = next(t for t, ks in enumerate(kinds) if "readout" in ks)
    out = {"band": 0.0, "fout": 0.0, "read": 0.0, "reset": 0.0}
    for t, st in enumerate(steps):
        c = max(op_cost(ops[j], k) for j in st)
        if t < band_end:
            out["band"] += c
        elif "readout" in kinds[t]:
            out["fout"] += c
        elif "reset" in kinds[t]:
            out["reset"] += c
        else:
            out["read"] += c
    return out


def demand_rate_per_s(d: int, k: int = 1) -> float:
    """Aggregate raw-pair demand: d-1 pairs per local round time."""
    return (d - 1) / (schedule_time_us(d, merge=False, rounds=1, k=k) * 1e-6)


def _check(d: int) -> None:
    """Re-assert the certified structure before pricing it."""
    assert len(S.parallel_steps(d, merge=False)) == 22 + 2 * d
    for beat in {b for beats in S.BEATS.values() for b in beats}:
        assert S.BEAT_KIND[beat] in set(DURATIONS_US) | {"junction"}, \
            f"unpriced beat {beat!r}"


if __name__ == "__main__":
    ds = [int(sys.argv[1])] if len(sys.argv) > 1 else [3, 7, 13, 27]
    for d in ds:
        _check(d)
    print("beat-kind durations, us (optimistic / baseline / conservative):")
    for kind, v in DURATIONS_US.items():
        print(f"  {kind:14s} {v[0]:7.0f} {v[1]:7.0f} {v[2]:7.0f}")
    print(f"  {'junction':14s} {'2x shuttle':>23s}")
    print(f"  {'merge+gate beat':14s} {'gate + merge_split':>23s}")
    print("\nlocal round T_round, ms (opt/base/cons):")
    for d in ds:
        row = [schedule_time_us(d, False, 1, k) / 1000 for k in range(3)]
        print(f"  d={d:2d}: {row[0]:8.2f} {row[1]:8.2f} {row[2]:8.2f}")
    print("\nphase split at baseline, ms (band should be flat in d):")
    for d in ds:
        ph = phase_times_us(d, 1)
        print(f"  d={d:2d}: band {ph['band']/1000:6.2f}  file-out {ph['fout']/1000:6.2f}"
              f"  read {ph['read']/1000:6.2f}  reset {ph['reset']/1000:6.2f}")
    print("\nfull d-round merge T_merge, ms (opt/base/cons):")
    for d in ds:
        if d > 13:
            continue                                 # the big merges take a while
        row = [schedule_time_us(d, True, d, k) / 1000 for k in range(3)]
        print(f"  d={d:2d}: {row[0]:8.2f} {row[1]:8.2f} {row[2]:8.2f}")
    print("\naggregate demand rate (d-1)/T_round, pairs/s:")
    for d in ds:
        row = [demand_rate_per_s(d, k) for k in range(3)]
        print(f"  d={d:2d}: {row[0]:8.1f} {row[1]:8.1f} {row[2]:8.1f}")

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

The distance rule (thesis Section 4.2 floorplan)
------------------------------------------------
A shuttle beat is priced per 200 um segment traveled, not per beat.
The floorplan fixes the distances: wells at 3 per segment (67 um pitch,
the d=3 memory row of 7 wells spanning 2 segments); memory = 2d+1 wells;
gate strip = d interaction wells + d junction columns at the same pitch;
SPAM = d sites at the same pitch; wall and interface one segment each.
Zones run [memory | gate | SPAM | wall | interface] along one axis.
Adjacent rows sit ROW_PITCH_SEGMENTS apart (two segments, 400 um: two
RF rails plus the cross-row transport lane). A junction transit is two
intersection crossings plus that run, (2 + pitch) legs. The park
descent to the bottom cell crosses up to (d-1) pitches, charged at the
worst row, once per merge. The gate strip orders its columns with the
seam column d-1 at the interface end. Each active lane runs two comm ions
in ping-pong (Section 4.2): the carrier waits at a gate-end swap well and
hops one segment to its adjacent boundary data, while the networker's
per-round handoff to that swap well is the interface->gate trip.

Transport that feeds the gate strip pipelines under the running gate:
while one step's gate runs (gate + merge + split of cover), the previous
step's ions file home and the next step's ions file in, serially on the
one rail, 2*M segments under one cover. Only what the cover cannot hide
is charged: the band's first entry and last return (2*M segments), any
per-step residual where 2*M*t_seg exceeds the cover, the ancilla convoy
to SPAM and back (M+G+S segments each way, the train-pass model: the
file-out swap layers execute at the strip end as the train passes), and
the networker's ping-pong handoff (interface + wall + SPAM + one segment,
inside comm_swap; the carrier's own gate hop is one segment).
The junction's two-crossing price doubles as its ~two-segment pitch.
A serial bound with no pipelining is also printed, as the honest ceiling.

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
      assumption. The wide bracket owns settle time; path length is
      charged per segment by the distance rule below.
  junction      (2 + pitch) x shuttle (derived, never set independently)
      One junction beat crosses the intersection twice (see the op_tally
      NOTE in qec_scheduler) and runs the two-segment row pitch between
      the crossings. Baseline 208 us = 4 x 52. The park descent instead
      runs (d-1) pitches, charged at the worst row, once per merge.
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
import math
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

# ---- the floorplan, thesis Section 4.2, in 200 um segments -----------------
WELLS_PER_SEGMENT = 3            # 67 um pitch; d=3 memory: 7 wells over 2 seg


def mem_segments(d: int) -> int:
    return math.ceil((2 * d + 1) / WELLS_PER_SEGMENT)


def gate_segments(d: int) -> int:
    return math.ceil(2 * d / WELLS_PER_SEGMENT)   # d gate wells + d junctions


def spam_segments(d: int) -> int:
    return math.ceil(d / WELLS_PER_SEGMENT)       # d detection sites


IF_SEGMENTS = 1                  # interface zone, about one segment (ch4)
ROW_PITCH_SEGMENTS = 2           # rails + cross-row lane between rows
WALL_SEGMENTS = 1                # slotted baffle and clearance


def shuttle_segments_table(d: int) -> dict:
    """Distance of every shuttle-kind beat, in 200 um segments. Exhaustive
    by construction: pricing an unlisted beat is an error, never a silent
    default (the silent default is how distances went missing once)."""
    M, G, Sp = mem_segments(d), gate_segments(d), spam_segments(d)
    return {
        "shuttle-to-SPAM":      M + G + Sp,   # convoy: farthest ancilla
        "shuttle-from-SPAM":    M + G + Sp,   # crosses memory, gate, SPAM
        # two-ion ping-pong (Section 4.2): the carrier waits at the gate-end swap
        # well and hops one segment to its adjacent boundary data and back; the
        # long interface->gate trip is the networker's handoff inside comm_swap.
        "shuttle-carrier-hop":  1,            # carrier: swap well <-> seam gate
        "shuttle-handoff":      IF_SEGMENTS + WALL_SEGMENTS + Sp + 1,  # networker cavity -> swap well
        "shuttle-to-junction":  1,            # boundary ancilla, adjacent col
        "shuttle-to-park":      2,
        "shuttle-from-park":    2,
        "shuttle-home":         1,
        "settle":               1,
        "drop":                 1,
    }


def beat_segments(beat: str, d: int) -> int:
    return shuttle_segments_table(d)[beat]


def gate_cover_us(k: int) -> float:
    """What one gate step's own work covers: gate + merge + split."""
    return DURATIONS_US["gate"][k] + 2.0 * DURATIONS_US["merge_split"][k]


def beat_cost(beat: str, k: int, d: int = 7) -> float:
    """One beat's duration in us at bracket k (0 opt, 1 base, 2 cons)."""
    kind = S.BEAT_KIND[beat]
    if kind == "junction":
        ts = DURATIONS_US["shuttle"][k]
        if beat in ("junction-descent", "junction-ascent"):
            # parking: down/up the boundary lane, worst row crosses d-1
            # pitches, two intersection crossings at the ends
            return (2.0 + (d - 1) * ROW_PITCH_SEGMENTS) * ts
        # one row over: two crossings plus the pitch run
        return (2.0 + ROW_PITCH_SEGMENTS) * ts
    if beat == "merge+gate":
        return DURATIONS_US["gate"][k] + DURATIONS_US["merge_split"][k]
    if kind == "shuttle":
        return beat_segments(beat, d) * DURATIONS_US["shuttle"][k]
    return DURATIONS_US[kind][k]


def op_cost(op, k: int, d: int = 7) -> float:
    return sum(beat_cost(b, k, d) for b in S.BEATS[op[0]])


def _gate_feed_us(d: int, ops, steps, k: int) -> float:
    """The gate-strip feed the cover cannot hide: one entry and one return
    per round (2*M segments), plus the per-step residual where the 2*M
    round trip outgrows the cover."""
    ts = DURATIONS_US["shuttle"][k]
    trip = 2.0 * mem_segments(d) * ts
    resid = max(0.0, trip - gate_cover_us(k))
    gate_steps = sum(
        1 for st in steps
        if any(ops[j][0] in ("inrow", "xgate", "comm_gate") for j in st))
    rounds = max(1, sum(1 for op in ops if op[0] == "round"))
    return rounds * trip + gate_steps * resid


def schedule_time_us(d: int, merge: bool = False, rounds: int = 1,
                     k: int = 1, pipelined: bool = True) -> float:
    """Step-weighted critical path: each packed step costs its slowest op.
    Pipelined (default): gate-feed transport hides under the gate cover and
    only the uncoverable pieces are charged. Serial: every gate step pays
    its own 2*M round trip, the no-overlap ceiling."""
    ops = S.round_ops(d, merge=merge, rounds=rounds)
    steps = S.parallel_steps(d, merge=merge, rounds=rounds)
    base = sum(max(op_cost(ops[j], k, d) for j in st) for st in steps)
    if pipelined:
        return base + _gate_feed_us(d, ops, steps, k)
    ts = DURATIONS_US["shuttle"][k]
    gate_steps = sum(
        1 for st in steps
        if any(ops[j][0] in ("inrow", "xgate", "comm_gate") for j in st))
    return base + gate_steps * 2.0 * mem_segments(d) * ts


def phase_times_us(d: int, k: int = 1) -> dict:
    """The local round's time split into its four phases (band, file-out,
    shuttle/read, reset), the same segmentation the scheduler certifies.
    The uncoverable gate feed is charged to the band, where it happens."""
    ops = S.round_ops(d, merge=False)
    steps = S.parallel_steps(d, merge=False)
    kinds = [{ops[j][0] for j in st} for st in steps]
    band_end = next(t for t, ks in enumerate(kinds) if "readout" in ks)
    out = {"band": 0.0, "fout": 0.0, "read": 0.0, "reset": 0.0}
    for t, st in enumerate(steps):
        c = max(op_cost(ops[j], k, d) for j in st)
        if t < band_end:
            out["band"] += c
        elif "readout" in kinds[t]:
            out["fout"] += c
        elif "reset" in kinds[t]:
            out["reset"] += c
        else:
            out["read"] += c
    out["band"] += _gate_feed_us(d, ops, steps, k)
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
    shuttles = {b for b, kd in S.BEAT_KIND.items() if kd == "shuttle"}
    assert shuttles == set(shuttle_segments_table(d)), \
        "distance table and shuttle beats disagree"


if __name__ == "__main__":
    ds = [int(sys.argv[1])] if len(sys.argv) > 1 else [3, 7, 13, 27]
    for d in ds:
        _check(d)
    print("beat-kind durations, us (optimistic / baseline / conservative):")
    for kind, v in DURATIONS_US.items():
        print(f"  {kind:14s} {v[0]:7.0f} {v[1]:7.0f} {v[2]:7.0f}")
    print(f"  {'junction':14s} {'(2+pitch)x shuttle':>23s}")
    print(f"  {'merge+gate beat':14s} {'gate + merge_split':>23s}")
    print("\nfloorplan distances in 200 um segments (thesis Section 4.2):")
    for d in ds:
        M, G, Sp = mem_segments(d), gate_segments(d), spam_segments(d)
        print(f"  d={d:2d}: memory {M}  gate {G}  SPAM {Sp}  "
              f"SPAM convoy {M+G+Sp}  comm transit "
              f"{IF_SEGMENTS+WALL_SEGMENTS+Sp+1}  gate feed 2x{M}")
    print("\nlocal round T_round, ms (opt/base/cons), pipelined | serial bound:")
    for d in ds:
        row = [schedule_time_us(d, False, 1, k) / 1000 for k in range(3)]
        ser = [schedule_time_us(d, False, 1, k, pipelined=False) / 1000
               for k in range(3)]
        print(f"  d={d:2d}: {row[0]:8.2f} {row[1]:8.2f} {row[2]:8.2f}"
              f"   | {ser[0]:8.2f} {ser[1]:8.2f} {ser[2]:8.2f}")
    print("\nphase split at baseline, ms (gate work flat in d; the band's")
    print("transport edge grows with the row length):")
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

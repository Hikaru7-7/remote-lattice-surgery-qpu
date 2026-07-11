# remote-lattice-surgery-qpu

The software of a modular trapped-ion QPU designed for remote lattice
surgery. The scheduler builds one full round of surface-code error
correction, and a complete remote lattice-surgery merge across two modules,
as an ordered list of physical operations on a segmented ion trap. Then it
checks, at every code distance, that the schedule is legal. Ions cannot pass in their one-dimensional channels,
a well holds one ion except the pair merged for a gate, junctions are exclusive,
and every gate fires in a well, never on a junction.

A timing layer prices the certified schedule in seconds from published
demonstrations, and a requirement layer inverts the priced demand into the
efficiency the design's ion--photon interface must deliver.

This is the software behind Sections 4.4, 5.2, 5.3, 6.3, 7.2, and 7.3 of the
thesis *Trapped-Ion Multi-Computer Design for Remote Lattice Surgery* (Keio
University, 2026). Every count and every claim in those sections can be
reproduced here.

**What it is not.** The scheduler moves ions, not quantum states. It proves
the round can be scheduled without conflict, and the timing layer prices it.
Neither returns a logical error rate. A separate injection layer does measure
logical error, but on the abstract surface and repetition codes, not on the
physical ion schedule. Noise on the physical schedule itself is still future
work.

## Requirements

Python 3. The scheduler, timing, requirement, fidelity, supply, distillation,
and visualizer layers use only the standard library. Two numerical layers need
add-ons: `qec_inject.py` uses `numpy`, and `qec_inject_stim.py` /
`qec_inject_stim_hw.py` use `numpy`, `stim`, and `pymatching`. Install those
with `pip install numpy stim pymatching`.

## Run it

```
python3 qec_scheduler.py        # named checks, then every check at every odd d = 3..27 (expect 21 PASS, 22 at d = 3)
python3 qec_scheduler.py 7      # one-distance report: placement, seam, depth, tally
python3 qec_visualizer.py      # build the d=3 and d=5 HTML animations
python3 qec_visualizer.py 7    # build the d=7 HTML animations
python3 qec_visualizer.py all  # rebuild every odd d = 3..27 and print the sweep table
python3 qec_timing.py          # price the schedule in seconds: durations, T_round, T_merge, demand rate
python3 make_requirement.py    # invert the demand into the required interface efficiency
python3 make_fidelity.py       # the seam-grade floor, the memory line, the required visibility
python3 make_supply.py         # the MEMS supply side: geometry to the fiber, scenarios, verdict preview
python3 qec_distill.py         # certify the double-selection distillation round at every odd d
python3 qec_inject.py          # code-capacity injection: the space-like seam effect (numpy)
python3 qec_inject_stim.py     # circuit-level seam factor and its recovery (needs stim + pymatching)
python3 qec_inject_stim_hw.py  # circuit-level at the demonstrated hardware rates (needs stim + pymatching)
```

The default run ends by certifying every odd distance up to 27, so the claim
in the thesis is exactly what the command shows. One distance at a time is
`python3 qec_scheduler.py 27`.

## Files

| File | What it is |
|---|---|
| `qec_scheduler.py` | The source of truth. Builds the distance-d rotated surface code, places every ion on the chip, emits the schedule (`round_ops`), packs it into parallel time-steps (`parallel_steps`), and runs the structural checks. `op_tally` counts every physical beat so a duration model can turn the schedule into a round time. |
| `qec_timing.py` | Prices the certified schedule in seconds. Each packed step costs its slowest operation; per-beat durations are traced to published demonstrations and bracketed optimistic/baseline/conservative. Every round-time and demand-rate number in the thesis reprints from this file. |
| `make_requirement.py` | Inverts the priced demand into the required ion--photon interface efficiency, in the mean and the 99%-delivery forms. Chain brackets mirror Table 5.1 of the thesis; the timing comes live from `qec_timing.py`. Every requirement number in thesis Section 5.2 reprints from this file. |
| `make_fidelity.py` | The fidelity track of thesis 5.3: the seam-grade floor from the published seam tolerance, the memory line, and the required two-photon visibility. |
| `make_supply.py` | The supply side of thesis Chapter 6: cavity geometry to interface efficiency at the fiber, the three scenarios, and the verdict preview against the requirement. |
| `qec_distill.py` | The distillation track of thesis Section 7.3: certifies the double-selection round at every odd distance, and prices its pair cost and rate factor (3 raw pairs to 1, demand times 3.4). |
| `qec_inject.py` | Code-capacity error injection on the rotated surface code and the lattice-surgery seam, decoded by exact minimum-weight matching. Measures the space-like seam effect, about 1.5x the bulk rate. Uses `numpy`. |
| `qec_inject_stim.py` | Circuit-level error injection with Stim and PyMatching. Measures the merge's time-like logical error, the seam factor `1.8^((d+1)/2)` (about 11 at d = 7), and its recovery by one distance step. Uses `numpy`, `stim`, `pymatching`. |
| `qec_inject_stim_hw.py` | The same circuit-level measurement pinned to the demonstrated per-operation error rates. Uses `numpy`, `stim`, `pymatching`. |
| `CIRCUIT_LEVEL_RESULTS.md` | The captured output of the circuit-level runs, tabulated for cross-checking against thesis Sections 5.3 and 7.2. |
| `qec_visualizer.py` | Replays the scheduler's operation list, unchanged, on the full Chapter-4 cell geometry: the 2d+1 memory homes, the gate strip with its junction columns, gate wells, hold wells, and swap well, the SPAM sites, the wall, and the cavities, one cell per code row. The scheduler's `frame_errors` legality check runs on every frame, at build time and precomputed into each page's live line, so the visualizer reflects the scheduler and adds no rule of its own: well occupancy within capacity, no ion at rest on a junction column, and no reordering along a row without a shared-well crystal rotation. |
| `qec_round_sim_d{3,5,7}.html` | One local error-correction round at that distance, on the chip geometry. |
| `qec_merge_full_sim_d{3,5,7}.html` | Two rounds of the remote lattice-surgery merge, the seam walked by the comm ions. The full d-round merge is packed and certified in the scheduler. |
| `qec_merge_distill_sim_d{3,5,7}.html` | The same merge on the distilled lane, thesis Section 4.3.4: two halves ferried to the gate-end hold wells during the round, the spent survivor recycled at the read to catch the third, double selection at the boundary. Op order from `qec_distill.py`. |

## Viewing the animations

Watch them in the browser, nothing to download:

| | d = 3 | d = 5 | d = 7 |
|---|---|---|---|
| One local round | [round d3](https://hikaru7-7.github.io/remote-lattice-surgery-qpu/qec_round_sim_d3.html) | [round d5](https://hikaru7-7.github.io/remote-lattice-surgery-qpu/qec_round_sim_d5.html) | [round d7](https://hikaru7-7.github.io/remote-lattice-surgery-qpu/qec_round_sim_d7.html) |
| Full merge | [merge d3](https://hikaru7-7.github.io/remote-lattice-surgery-qpu/qec_merge_full_sim_d3.html) | [merge d5](https://hikaru7-7.github.io/remote-lattice-surgery-qpu/qec_merge_full_sim_d5.html) | [merge d7](https://hikaru7-7.github.io/remote-lattice-surgery-qpu/qec_merge_full_sim_d7.html) |
| Merge with distillation | [distilled d3](https://hikaru7-7.github.io/remote-lattice-surgery-qpu/qec_merge_distill_sim_d3.html) | [distilled d5](https://hikaru7-7.github.io/remote-lattice-surgery-qpu/qec_merge_distill_sim_d5.html) | [distilled d7](https://hikaru7-7.github.io/remote-lattice-surgery-qpu/qec_merge_distill_sim_d7.html) |

Or open any HTML file locally in a browser. No server needed. Use Prev / Next, the
slider, or Play. Every page draws the full cell geometry and every frame shows
a caption of the physical operation plus a live verifier line that turns red on
any violation: a well over its capacity, an ion at rest on a junction column,
or two ions changing order without a shared-well rotation, an ion resting on
a junction column, two ions sharing one junction point, a mid-transit change
of junction column, an ion entering a row anywhere but at a junction
column's channel, a gate firing on anything but an isolated pair in a gate
well, a read away from a SPAM site or the swap well, a communication ion in
the memory zone, or a code ion past the optical wall. It never fires. Every
move is a continuous path: along the row through explicit pairwise rotations,
vertical only at a junction column. Occupancy law: one ion rests per well,
gates fire on isolated pairs, parking in the gate strip reaches three, one
passing transit may briefly make four (the transport crystals of Pino et al.
2021), and the swap well never exceeds its resident pair. The
verdict is the scheduler's own, computed per frame and shipped to the page, so
the reader can see it pass rather than take it on trust. The capacities are the design's own:
one ion rests per well, gates and swaps merge pairs, a gate well carries at
most the four-ion crystal of thesis Table 4.1, and a swap well carries three
where the distilled lane's survivor waits beside the ping-pong pair. The
frame census across all hosted pages never exceeds them.

## The full sweep

The thesis claim covers every odd distance from 3 to 27, and the hosted pages
above are samples of it. One command re-verifies all three families at every
swept distance (a slot-level check, no pages written, a few minutes) and
reprints the table below; `python3 qec_visualizer.py 9` writes any one
distance's pages locally:

```
python3 qec_visualizer.py all
```

Every verifier rule also carries a failing test. `python3
qec_visualizer.py selftest` plants one violation of each rule into a clean
run, a teleport, a free pass, an over-capacity well, a rest on a junction
column, a shared junction point, and a mid-transit column switch, and
requires the verifier to refuse each one. A rule without a failing test is
a blind spot waiting to happen.

Each row runs the scheduler's named checks and its per-frame legality pass
(`frame_errors`) over the visualizer's frames at one distance. Local round steps are exactly `22 + 2d` at every
distance, and no frame at any distance produces a finding: capacities hold,
junction columns stay clear, and nothing passes without a rotation. The d = 3
row runs one extra check, the hand-verified reference comparison.

| d | scheduler checks | local round steps | frames | findings | 2-round merge steps | frames | findings | distilled merge frames | findings |
|--:|:----------------:|------------------:|-------:|---------:|--------------------:|-------:|---------:|-----------------------:|---------:|
| 3 | 22 PASS | 28 | 270 | 0 | 52 | 613 | 0 | 1010 | 0 |
| 5 | 21 PASS | 32 | 1063 | 0 | 64 | 2283 | 0 | 3170 | 0 |
| 7 | 21 PASS | 36 | 2549 | 0 | 72 | 5495 | 0 | 6968 | 0 |
| 9 | 21 PASS | 40 | 4893 | 0 | 80 | 10729 | 0 | 12884 | 0 |
| 11 | 21 PASS | 44 | 8218 | 0 | 88 | 17967 | 0 | 20900 | 0 |
| 13 | 21 PASS | 48 | 12586 | 0 | 96 | 27735 | 0 | 31542 | 0 |
| 15 | 21 PASS | 52 | 18188 | 0 | 104 | 40027 | 0 | 44804 | 0 |
| 17 | 21 PASS | 56 | 25168 | 0 | 112 | 55565 | 0 | 61408 | 0 |
| 19 | 21 PASS | 60 | 33672 | 0 | 120 | 74021 | 0 | 81026 | 0 |
| 21 | 21 PASS | 64 | 43844 | 0 | 128 | 96589 | 0 | 104852 | 0 |
| 23 | 21 PASS | 68 | 55828 | 0 | 136 | 122417 | 0 | 132034 | 0 |
| 25 | 21 PASS | 72 | 69768 | 0 | 144 | 153245 | 0 | 164312 | 0 |
| 27 | 21 PASS | 76 | 85808 | 0 | 152 | 187613 | 0 | 200226 | 0 |

The full sweep, every odd distance from 3 to 27, all three families, zero
findings at every row: some 1.5 million frames of schedule under the eleven
rules. The command reprints this table in a few minutes of CPU.

## Headline numbers to reproduce

- Local round depth is exactly `22 + 2d` parallel time-steps at every distance.
  A flat 19-step connection band, a readout tail of `d + 3`, then `d` swap
  layers back to the rest state. This is an asserted check, not a printout.
- A cross-row ancilla swaps toward its junction only when a data ion blocks
  the way; the junction-adjacent half of the crossings lift directly, and
  those free lifts overlap the in-row gates. That is why the band is 19.
- Every schedule ends in the placement it started from, every swap replayed
  and undone. This is its own check, so rounds and merges chain.
- The full d-round merge packs into 78, 160, and 252 time-steps at d = 3, 5, 7,
  and every one of its d rounds takes exactly the same number of steps.
- Resources are closed forms. `d^2` data ions, `d^2 - 1` ancillas, `d` comm
  lanes, `d - 1` Bell pairs per merge round, `4d(d - 1)` two-qubit gates per
  round, `floor((d-1)/2)` park wells.

## Design contract

The scheduler owns the schedule. What happens, in what order, with what motion,
is decided in `qec_scheduler.py` and certified there. The visualizer expands that op list to per-frame ion positions and hands each
frame to `qec_scheduler.frame_errors`; it draws the schedule and reflects that
verdict, adding no rule of its own. It still rejects any operation it has no
renderer for.

## The circuit-level layer

The scheduler is certified and priced. A separate layer injects noise and
measures logical error, so the seam tolerance of remote lattice surgery is
measured here rather than cited from the literature.

- `qec_inject.py` runs code-capacity noise on the rotated surface code and on
  the lattice-surgery seam, decodes with exact minimum-weight matching, and
  measures the space-like seam effect, about 1.5x the bulk rate and roughly
  flat in distance.
- `qec_inject_stim.py` and `qec_inject_stim_hw.py` run full circuit-level noise
  with Stim and PyMatching. They measure the merge's time-like logical error,
  the seam factor `1.8^((d+1)/2)` (about 11 at d = 7), and its recovery by one
  distance step. The `_hw` variant pins every noise knob to the demonstrated
  per-operation rate.
- `qec_distill.py` certifies the double-selection distillation round at every
  odd distance.

These run on the abstract code circuits, not on the physical ion schedule;
injecting noise on the physical schedule itself remains future work. The
captured results are in `CIRCUIT_LEVEL_RESULTS.md`.

## Cite

```bibtex
@misc{Yokomori2026QECScheduler,
  author       = {Yokomori, Hikaru},
  title        = {remote-lattice-surgery-qpu: certified scheduling, timing, and interface requirements for a modular trapped-ion {QPU}},
  year         = {2026},
  howpublished = {\url{https://github.com/Hikaru7-7/remote-lattice-surgery-qpu}},
  note         = {Open-source software}
}
```

## License

MIT. See `LICENSE`.

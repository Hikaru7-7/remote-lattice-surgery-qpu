# Circuit-level error-injection results (Stim + PyMatching)

Captured record of the circuit-level seam measurements for thesis Section 5.3 / 7.2.
These were run with Stim 1.15.0 and PyMatching 2.4.0 on a machine where those wheels
install (the build sandbox is Linux ARM64 + Python 3.10, which Stim ships no wheel for).
Reproduce with `python3 qec_inject_stim.py` and `python3 qec_inject_stim_hw.py`.

Noise model: circuit-level (every gate, reset, and measurement can fail). Decoding is
minimum-weight matching off Stim's detector error model. The merge's time-like logical
(the joint parity held over d rounds) is a repetition code in the round direction; the
seam carries the eight charged operations of the merge, modelled as 1.8x the bulk rate.

## 1. The seam factor vs distance (uniform p = 0.007, 3,000,000 shots)

Same code at 1.8x vs 1x seam noise; the ratio is the measured factor.

| d | 1x seam | 1.8x seam | measured factor | analytic 1.8^((d+1)/2) | 1x fails |
|---|---------|-----------|-----------------|------------------------|----------|
| 3 | 3.68e-3 | 1.12e-2   | 3.05            | 3.2                    | 11027    |
| 5 | 5.06e-4 | 2.90e-3   | 5.74            | 5.8                    | 1518     |
| 7 | 7.0e-5  | 7.49e-4   | 10.70           | 10.5                   | 210      |
| 9 | 1.0e-5  | 1.79e-4   | 17.35           | 18.9                   | 31       |

The measured factor tracks 1.8^((d+1)/2); at d=7 it is ~10.7, matching the analytic 10.5.

## 2. Recovery: does one distance step undo the factor? (5,000,000 shots)

Loaded (1.8x) merge time-like; per-step suppression = loaded(7)/loaded(9).

| p     | loaded(7) | loaded(9) | 1-step suppression | supp x p |
|-------|-----------|-----------|--------------------|----------|
| 0.004 | 7.54e-5   | 1.48e-5   | 5.1x               | 0.020    |
| 0.005 | 1.85e-4   | 3.44e-5   | 5.4x               | 0.027    |
| 0.006 | 4.03e-4   | 8.36e-5   | 4.8x               | 0.029    |
| 0.007 | 7.53e-4   | 1.81e-4   | 4.2x               | 0.029    |

`supp x p` is roughly constant, so suppression ~ threshold / (1.8 p), threshold ~ 0.05.
The suppression exceeds the ~10.5 factor below p ~ 0.003. The operating point (p ~ 1e-3)
sits below that, so one distance step recovers the factor at the operating error rate; it
only looks like two steps at the near-threshold rates where failures are directly countable.

## 3. At the demonstrated per-operation rates (3,000,000 shots)

Gate 1.8e-3 (Loeschnauer 2024, benchmarked average), SPAM 4.8e-4 (Ransford 2025 / Helios),
idle 2.8e-4 (memory line). Seam = 1.8x these.

Bulk rotated-surface-code memory logical error:

| d | p_L      | fails / 3,000,000 |
|---|----------|-------------------|
| 3 | 1.00e-3  | 3003              |
| 5 | 2.81e-4  | 843               |
| 7 | 5.67e-5  | 170               |
| 9 | 1.10e-5  | 33                |

Merge time-like logical error (unloaded = bulk rate, loaded = 1.8x):

| d | unloaded | loaded  | seam factor | analytic | loaded fails |
|---|----------|---------|-------------|----------|--------------|
| 3 | 6.93e-5  | 2.32e-4 | 3.35        | 3.2      | 696          |
| 5 | 2.0e-6   | 1.03e-5 | 5.17        | 5.8      | 31           |
| 7 | <1/3e6   | 3.3e-7  | (too rare)  | 10.5     | 1            |

At the demonstrated hardware the seam factor is confirmed at d=3 and d=5 (3.35, 5.17),
and the merge's time-like logical error at d=7 is below 1e-6 -- the seam, even at 1.8x the
bulk rate, is deep below threshold. The bulk memory logical error at d=7 is 5.7e-5.

#!/usr/bin/env python3
"""qec_inject_stim.py -- CIRCUIT-LEVEL seam measurement, thesis Section 7.2 / 5.3.

RUN THIS ON A MACHINE WHERE STIM INSTALLS (your Mac):
    pip3 install stim pymatching
    python3 qec_inject_stim.py
then paste the output back. (The build sandbox is Linux ARM64 + Python 3.10, a
combination Stim ships no wheel for, so this one script runs on your side.)

What it measures, at CIRCUIT LEVEL -- every gate, reset, and measurement can fail,
the honest noise model the code-capacity qec_inject.py could not reach:

  1. The bulk rotated-surface-code memory logical error, the reference scale.
  2. The merge's TIME-LIKE logical error -- the joint ZZ parity, held by measuring the
     seam checks every round for d rounds, is a repetition code in the round direction.
     Loading that seam at 1.8x the bulk rate (the eight charged operations of Section 5.3)
     and taking the ratio to the unloaded seam, SAME code, isolates the seam factor. Deep
     below threshold that ratio is 1.8^((d+1)/2), which is ~8 at d=7 -- the analytical
     claim, now measured rather than cited.

Decoding is minimum-weight matching (PyMatching) off Stim's detector error model, so no
decoder tuning enters the number.
"""
import numpy as np
import stim
import pymatching


def _tee_stdout(prefix):
    """Send everything printed to results/<prefix>_<timestamp>.txt as well as the terminal,
    so each run is captured to a file. Returns the path."""
    import sys, os, datetime
    os.makedirs("results", exist_ok=True)
    path = os.path.join("results", f"{prefix}_{datetime.datetime.now():%Y%m%d_%H%M%S}.txt")
    fh, real = open(path, "w"), sys.stdout

    class _Tee:
        def write(self, s):
            real.write(s); fh.write(s)

        def flush(self):
            real.flush(); fh.flush()
    sys.stdout = _Tee()
    return path


def surface_memory(d, p, rounds=None):
    return stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        distance=d, rounds=rounds or d,
        after_clifford_depolarization=p,
        before_round_data_depolarization=p,
        before_measure_flip_probability=p,
        after_reset_flip_probability=p)


def seam_timelike(d, p, rounds=None):
    """The merge holds the joint parity by measuring the seam checks every round; the
    time-like logical is that parity, protected as a repetition code of distance = rounds."""
    return stim.Circuit.generated(
        "repetition_code:memory",
        distance=d, rounds=rounds or d,
        after_clifford_depolarization=p,
        before_round_data_depolarization=p,
        before_measure_flip_probability=p,
        after_reset_flip_probability=p)


def logical_error_rate(circuit, shots):
    dem = circuit.detector_error_model(decompose_errors=True,
                                       approximate_disjoint_errors=True)
    matcher = pymatching.Matching.from_detector_error_model(dem)
    sampler = circuit.compile_detector_sampler()
    dets, obs = sampler.sample(shots, separate_observables=True)
    preds = matcher.decode_batch(dets)
    return np.count_nonzero(np.any(preds != obs, axis=1)) / shots


if __name__ == "__main__":
    _saved = _tee_stdout("qec_inject_stim")
    print(f"qec_inject_stim.py  (stim {stim.__version__} | pymatching {pymatching.__version__})\n")
    SHOTS = 300_000
    SEAM = 1.8   # seam / bulk error ratio, the eight charged operations (Section 5.3)
    DS = (3, 5, 7, 9)

    print("1) bulk rotated-surface-code memory, circuit-level logical error:")
    print(f"   {'p':>7} | " + " ".join(f"d={d:<2d}    " for d in DS))
    for p in (0.001, 0.002, 0.003, 0.005):
        row = [logical_error_rate(surface_memory(d, p), SHOTS) for d in DS]
        print(f"   {p:>7.3f} | " + " ".join(f"{x:9.6f}" for x in row))

    print("\n2) merge time-like logical error (seam checks over d rounds), seam at 1.8x:")
    print(f"   {'p':>7} | " + " ".join(f"d={d:<2d}    " for d in DS))
    for p in (0.001, 0.002, 0.003, 0.005):
        row = [logical_error_rate(seam_timelike(d, SEAM * p), SHOTS) for d in DS]
        print(f"   {p:>7.3f} | " + " ".join(f"{x:9.6f}" for x in row))

    print("\n3) THE SEAM FACTOR on the time-like logical (same code, 1.8x vs 1x noise):")
    print("   the measured version of the analytical 'near 8' at d=7. High statistics.")
    SHOTS_HI = 3_000_000
    p = 0.007
    print(f"   at p={p}, {SHOTS_HI:,} shots (fails column shows the statistics):")
    print(f"     d |   1x seam |  1.8x seam | factor | analytic 1.8^((d+1)/2) | 1x fails")
    for d in DS:
        base = logical_error_rate(seam_timelike(d, p), SHOTS_HI)
        load = logical_error_rate(seam_timelike(d, SEAM * p), SHOTS_HI)
        fac = load / base if base > 0 else float("nan")
        print(f"     {d:>2} | {base:9.6f} | {load:10.6f} | {fac:6.2f} | "
              f"{1.8 ** ((d + 1) / 2):20.1f} | {int(round(base * SHOTS_HI)):8d}")

    print("\n4) RECOVERY: does one distance step undo the seam factor at d=7?")
    print("   loaded = 1.8x seam. per-step suppression = loaded(7)/loaded(9);")
    print("   one step recovers the ~12x factor once this suppression exceeds ~12.")
    SHOTS_R = 5_000_000
    print(f"     p     |  loaded(7)   loaded(9) | 1-step supp | supp x p")
    for p in (0.004, 0.005, 0.006, 0.007):
        ld7 = logical_error_rate(seam_timelike(7, SEAM * p), SHOTS_R)
        ld9 = logical_error_rate(seam_timelike(9, SEAM * p), SHOTS_R)
        supp = ld7 / ld9 if ld9 > 0 else float("nan")
        print(f"   {p:.3f}  | {ld7:.3e}   {ld9:.3e} | {supp:7.1f}x  | {supp * p:.3f}")
    print("   (supp x p roughly constant => suppression = C/p; one step recovers where")
    print("    C/p >= 12, i.e. below p ~ C/12 -- compare that to the operating p~1e-3.)")

    print(f"\n[output saved to {_saved}]")
    print("--- paste everything above this line back ---")

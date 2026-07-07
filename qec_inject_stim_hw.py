#!/usr/bin/env python3
"""qec_inject_stim_hw.py -- circuit-level seam measurement at the DEMONSTRATED per-operation
error rates, thesis Section 7.2 / 5.3. Companion to qec_inject_stim.py, which scans a generic p.

RUN ON A MACHINE WHERE STIM INSTALLS (your Mac):
    pip3 install stim pymatching
    python3 qec_inject_stim_hw.py
then paste the output back.

Instead of a uniform depolarizing p, every Stim noise knob is set to the rate the thesis
already cites for that operation, so the merge's logical error and the seam factor reprint
from the demonstrations rather than from a chosen number:

  * two-qubit gate  1.8e-3   benchmarked-average microwave-gradient gate, Loeschnauer 2024
                             (best Bell 3e-4; the reference takes the average, Section 4.2)
  * SPAM (measure + reset)  4.8e-4   combined state-prep-and-measurement, Ransford 2025 (Helios)
  * data idle / round       2.8e-4   the memory line T_round/T2 at the baseline (make_fidelity.py)

The seam carries the eight charged operations of the merge (Section 5.3), modelled as 1.8x
these rates. Decoding is minimum-weight matching (PyMatching) off Stim's detector error model.
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


P_GATE = 1.8e-3   # Loeschnauer 2024, benchmarked-average two-qubit microwave-gradient gate
P_SPAM = 4.8e-4   # Ransford 2025 (Helios), combined state-preparation-and-measurement
P_IDLE = 2.8e-4   # memory line T_round/T2 at the baseline (make_fidelity.py)
SEAM = 1.8        # the eight charged operations of the seam, as a ratio to the bulk (Section 5.3)


def _gen(kind, d, scale):
    return stim.Circuit.generated(
        kind, distance=d, rounds=d,
        after_clifford_depolarization=P_GATE * scale,
        before_measure_flip_probability=P_SPAM * scale,
        after_reset_flip_probability=P_SPAM * scale,
        before_round_data_depolarization=P_IDLE * scale)


def surface_hw(d, scale=1.0):
    return _gen("surface_code:rotated_memory_z", d, scale)


def seam_hw(d, scale=1.0):
    """The merge's time-like logical: the joint parity held over d rounds, a repetition code
    in the round direction, at the demonstrated rates times scale."""
    return _gen("repetition_code:memory", d, scale)


def logical_error_rate(circuit, shots):
    dem = circuit.detector_error_model(decompose_errors=True,
                                       approximate_disjoint_errors=True)
    matcher = pymatching.Matching.from_detector_error_model(dem)
    dets, obs = circuit.compile_detector_sampler().sample(shots, separate_observables=True)
    preds = matcher.decode_batch(dets)
    return np.count_nonzero(np.any(preds != obs, axis=1)) / shots


if __name__ == "__main__":
    _saved = _tee_stdout("qec_inject_stim_hw")
    print(f"qec_inject_stim_hw.py  (stim {stim.__version__} | pymatching {pymatching.__version__})")
    print("circuit-level noise at the DEMONSTRATED per-operation rates:")
    print(f"  gate {P_GATE:.1e} (Loeschnauer 2024)  SPAM {P_SPAM:.1e} (Ransford 2025)  "
          f"idle {P_IDLE:.1e} (memory line)")
    print(f"  seam = 1.8x these rates (eight charged operations, Section 5.3)\n")
    SHOTS = 3_000_000
    DS = (3, 5, 7, 9)

    print("1) bulk rotated-surface-code memory logical error, at the demonstrated rates:")
    for d in DS:
        b = logical_error_rate(surface_hw(d, 1.0), SHOTS)
        print(f"     d={d}:  p_L = {b:.3e}   ({int(round(b*SHOTS))} fails / {SHOTS:,})")

    print("\n2) merge time-like logical error at the demonstrated rates:")
    print("   (seam checks over d rounds; unloaded = bulk rate, loaded = 1.8x for the 8 ops)")
    print(f"     d |  unloaded   |   loaded    | seam factor | analytic 1.8^((d+1)/2) | loaded fails")
    for d in (3, 5, 7):
        un = logical_error_rate(seam_hw(d, 1.0), SHOTS)
        ld = logical_error_rate(seam_hw(d, SEAM), SHOTS)
        fac = ld / un if un > 0 else float("nan")
        print(f"     {d} | {un:.3e}  | {ld:.3e}  | {fac:10.2f} | "
              f"{1.8 ** ((d + 1) / 2):20.1f} | {int(round(ld*SHOTS)):8d}")

    print("\n   At the demonstrated bulk rate (1.8e-3, below the ~0.003 one-step-recovery")
    print("   crossover found in qec_inject_stim.py), one distance step recovers the factor.")
    print(f"\n[output saved to {_saved}]")
    print("--- paste everything above this line back ---")

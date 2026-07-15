#!/usr/bin/env python3
"""make_requirement.py -- the interface inversion of thesis 5.2 (rate side).

Success per attempt with symmetric arms:  p = (1/2) (eta_int * eta_chain)^2.
A lane fires N = T_round / T_attempt attempts inside its round window (the
pipelined herald shadow of Section 4.4.4). Two delivery forms are computed.

  mean form   p N >= 1        expect one pair per window
  0.99 form   P_del >= 0.99   with P_del = 1 - (1-p)^N, i.e. p N >= ln(100)

Both are inverted for the required interface efficiency eta_int, the
probability that one attempt leaves a photon in the lane's fiber. The chain
brackets mirror Table 5.1 of the thesis (sources there); the timing comes
live from qec_timing.py in this repository. Run:  python3 make_requirement.py
"""
import math
import sys
import qec_timing as T

BRACKETS = ("optimistic", "baseline", "conservative")

# chain elements, (opt, base, cons); sources in thesis Table 5.1
FIBER    = (0.99, 0.94, 0.71)       # 1 / 5 / 30 m at <= 50 dB/km, 493 nm
OPTICS   = (0.95, 0.90, 0.80)       # itemized allowance
DETECTOR = (0.80, 0.71, 0.60)       # SNSPD band edge / APD spec / APD low edge
T_ATT_US = (0.35, 1.0, 10.0)        # attempt cycle, us
BSA      = 0.5                      # exact heralded fraction, not a loss

def eta_chain(k):
    return FIBER[k] * OPTICS[k] * DETECTOR[k]

def required_eta_int(d, k, pn_target):
    """eta_int such that p*N = pn_target in bracket k at distance d."""
    t_window_us = T.merge_window_us(d, k)     # one amortized merge round
    n_attempts = t_window_us / T_ATT_US[k]
    p_req = pn_target / n_attempts
    eta_arm = math.sqrt(p_req / BSA)              # p = BSA * eta_arm^2
    return eta_arm / eta_chain(k), p_req, n_attempts

if __name__ == "__main__":
    d = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    # sanity: chain products match the thesis table to two digits
    for k, want in enumerate((0.75, 0.60, 0.34)):
        assert abs(eta_chain(k) - want) < 0.005, (k, eta_chain(k))
    print(f"chain products eta_chain: "
          + " ".join(f"{eta_chain(k):.3f}" for k in range(3)))
    ln100 = math.log(100.0)
    print(f"\nrequired interface efficiency eta_int at d={d}:")
    print(f"  {'bracket':14s} {'N/window':>9s} {'p_req(mean)':>12s} "
          f"{'eta_int(mean)':>14s} {'eta_int(99%)':>13s}")
    for k in range(3):
        e1, p1, n = required_eta_int(d, k, 1.0)
        e2, _, _ = required_eta_int(d, k, ln100)
        print(f"  {BRACKETS[k]:14s} {n:9.0f} {p1:12.3e} {100*e1:13.2f}% {100*e2:12.2f}%")

    # ---- 5.2 (link level): the per-attempt demand, both forms ----------
    # This is the demand BEFORE the chain is split off, the closing numbers
    # of thesis Section 5.2 and the per-attempt row of Table 5.4.
    print("\nper-attempt success demand p at the link (no chain split):")
    WANT_P99 = (2.42e-4, 3.30e-4, 9.56e-4)       # thesis quotes at the merge window
    for k in range(3):
        _, p_mean, n = required_eta_int(d, k, 1.0)
        _, p_99, _ = required_eta_int(d, k, ln100)
        print(f"  {BRACKETS[k]:14s} N={n:6.0f}  mean {p_mean:.2e}  "
              f"99% form {p_99:.2e}")
        assert abs(p_99 - WANT_P99[k]) / WANT_P99[k] < 0.05, (k, p_99)
    assert abs(required_eta_int(d, 1, 1.0)[1] - 7.17e-5) / 7.17e-5 < 0.02

    print(f"\ndelivery probability at exactly pN=1: {1 - math.exp(-1):.3f}")
    print(f"expected empty windows per merge, d*d(1-P): "
          f"{d*d*0.01:.2f} at P=0.99, {d*d*math.exp(-1):.1f} at pN=1")

    # ---- 5.4: which delivery form binds -------------------------------
    # A round fires when all d active lanes hold a pair (the seam has d checks:
    # d-1 weight-4 + 1 weight-2, one Bell pair each). Lanes that delivered keep
    # theirs, so the round waits E[W] windows, the mean of the maximum of d
    # geometric variables with per-window success P.
    def expected_windows(P, lanes):
        return sum(1.0 - (1.0 - (1.0 - P) ** k) ** lanes for k in range(200))
    for P, tag in ((1 - math.exp(-1), "mean form"), (0.99, "99% form")):
        w = expected_windows(P, d)
        print(f"  {tag:9s} P={P:.3f}: E[windows per round] = {w:.2f} "
              f"-> merge stretch x{w:.2f}")
    assert 3.0 < expected_windows(1 - math.exp(-1), 7) < 3.2   # d=7: 7 active lanes
    assert expected_windows(0.99, 7) < 1.08

    # ---- ch7 (discussion): the distilled operating point ---------------
    # One round of double selection consumes N_RAW_IN raw pairs and succeeds
    # when both checks pass, roughly (1-eps)^2 at raw error eps. The rate
    # demand scales by N_RAW_IN/P_ds and eta_int by its square root. N_RAW_IN
    # imported from qec_distill (the certified circuit) so the two cannot drift.
    import qec_distill as QD
    eps_raw = 0.06                       # baseline raw pair error (ch6)
    p_ds = (1 - eps_raw) ** 2
    factor = QD.N_RAW_IN / p_ds
    eta99_base = 100 * required_eta_int(d, 1, ln100)[0]
    print(f"\ndistilled point at baseline: P_ds ~ {p_ds:.2f}, "
          f"demand x{factor:.1f}, eta_int x{math.sqrt(factor):.2f} "
          f"-> {eta99_base*math.sqrt(factor):.1f}% (99% form)")
    assert 7.5 < eta99_base * math.sqrt(factor) < 9.5

    # ---- ch7 (discussion): the simple-distillation floor ---------------
    # Simple protocols floor near 10 p_loc (Campbell 2007; Krastanov 2019).
    # With the charged operations (make_fidelity.COMM_OPS), the seam ratio
    # against the bulk is (10 + COMM_OPS)/tol; COMM_OPS imported so the two
    # layers cannot drift (they used to hardcode 6/16 independently).
    import make_fidelity as MF
    seam_total = MF.MULTIPLIER + MF.COMM_OPS         # 10 + 8 = 18
    # time-like factor r^((d+1)/2), the thesis 5.3/7.2 formula; confirmed at
    # circuit level in CIRCUIT_LEVEL_RESULTS.md (measured 10.7 at d=7 vs 10.5)
    for tol in (10.0, 14.0):
        r = seam_total / tol
        print(f"  seam/bulk ratio after distillation at {tol:.0f}x: "
              f"{r:.2f} -> time-like factor {r**((d+1)/2):.1f} at d={d}")
    assert abs((seam_total/10.0)**((d+1)/2) - 10.5) < 0.1   # "about 11" (ch7)

    print(f"free-space comparator (demonstrated 2.3% per arm incl. 71% APD): "
          f"eta_int equiv {0.023/0.71*100:.1f}%")

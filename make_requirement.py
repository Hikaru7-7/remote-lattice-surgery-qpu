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
    t_round_us = T.schedule_time_us(d, merge=False, rounds=1, k=k)
    n_attempts = t_round_us / T_ATT_US[k]
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
    print(f"\ndelivery probability at exactly pN=1: {1 - math.exp(-1):.3f}")
    print(f"expected empty windows per merge, d(d-1)(1-P): "
          f"{d*(d-1)*0.01:.2f} at P=0.99, {d*(d-1)*math.exp(-1):.1f} at pN=1")

    # ---- 5.4: which delivery form binds -------------------------------
    # A round fires when all d-1 active lanes hold a pair. Lanes that
    # delivered keep theirs, so the round waits E[W] windows, the mean of
    # the maximum of d-1 geometric variables with per-window success P.
    def expected_windows(P, lanes):
        return sum(1.0 - (1.0 - (1.0 - P) ** k) ** lanes for k in range(200))
    for P, tag in ((1 - math.exp(-1), "mean form"), (0.99, "99% form")):
        w = expected_windows(P, d - 1)
        print(f"  {tag:9s} P={P:.3f}: E[windows per round] = {w:.2f} "
              f"-> merge stretch x{w:.2f}")
    assert 2.8 < expected_windows(1 - math.exp(-1), 6) < 3.1
    assert expected_windows(0.99, 6) < 1.07

    # ---- 5.4: the distilled operating point ----------------------------
    # One round of double selection consumes three raw pairs and succeeds
    # when both checks pass, roughly (1-eps)^2 at raw error eps. The rate
    # demand scales by 3/P_ds and eta_int by its square root.
    eps_raw = 0.06                       # baseline raw pair error (ch6)
    p_ds = (1 - eps_raw) ** 2
    factor = 3.0 / p_ds
    eta99_base = 100 * required_eta_int(d, 1, ln100)[0]
    print(f"\ndistilled point at baseline: P_ds ~ {p_ds:.2f}, "
          f"demand x{factor:.1f}, eta_int x{math.sqrt(factor):.2f} "
          f"-> {eta99_base*math.sqrt(factor):.1f}% (99% form)")
    assert 8.0 < eta99_base * math.sqrt(factor) < 9.5

    # ---- 5.4: the simple-distillation floor through the seam ----------
    # Simple protocols floor near 10 p_loc (Campbell 2007; Krastanov 2019).
    # With the six charged operations, the seam ratio against the bulk is
    # (10+6)/10 = 1.6 at the 10x tolerance, (10+6)/14 at the matched cut.
    for tol in (10.0, 14.0):
        r = 16.0 / tol
        print(f"  seam/bulk ratio after distillation at {tol:.0f}x: "
              f"{r:.2f} -> logical factor {r**(d/2):.1f} at d={d}")

    print(f"free-space comparator (demonstrated 2.3% per arm incl. 71% APD): "
          f"eta_int equiv {0.023/0.71*100:.1f}%")

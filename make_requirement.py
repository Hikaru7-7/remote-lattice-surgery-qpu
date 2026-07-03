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
    print(f"free-space comparator (demonstrated 2.3% per arm incl. 71% APD): "
          f"eta_int equiv {0.023/0.71*100:.1f}%")

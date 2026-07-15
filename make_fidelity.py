#!/usr/bin/env python3
"""make_fidelity.py -- the fidelity track of thesis Section 5.3.

The seam-grade floor. Ramette et al. (npj QI 10:58, 2024) show the merged
patch fails as two decoupled systems, bulk and seam, with thresholds near
1% per local gate and 10% per Bell pair. Holding the seam to the bulk's
standard gives eps_eff <= 10 p_loc. The seam check spends local work too:
two copy gates and one read per module each round, plus the ping-pong
reorder that hands the heralded half from the networking Ba+ to the
carrier Ba+ (Section 4.2, the two-comm-ion lane), one reorder per module.
That is four local operations per module, eight per check, charged in
full, so the pair's share is

    eps  <=  (10 - 8) p_loc  =  2 p_loc        (the seam-grade floor)

The memory line. The stored half of a pipelined pair waits at most one
round. Charged in full as T_round / T2 (upper bound, no decay-shape
argument), with T_round live from qec_timing.py and T2 anchored at the
clock-qubit demonstrations cited in the thesis (50 s Ca-43 with nothing,
10 min Yb-171 with dynamical decoupling).

Run:  python3 make_fidelity.py
"""
import sys
import qec_timing as T

D = 7
MULTIPLIER = 10          # Ramette Eq. (7) working number (14 = matched cut)
COMM_OPS = 8             # per module per check: 2 copy gates + 1 read + 1
                         # ping-pong reorder (2 comm ions/lane, Section 4.2)
P_LOC = {"family floor (Loeschnauer 2024)": 1e-3,
         "demonstrated 99.7% (Harty 2016)": 3e-3}
T2_S = {"Ca-43 clock, no decoupling (Harty 2014)": 50.0,
        "Yb-171 clock, decoupled (Wang 2017)": 600.0}
BRACKETS = ("optimistic", "baseline", "conservative")

if __name__ == "__main__":
    d = int(sys.argv[1]) if len(sys.argv) > 1 else D
    share = MULTIPLIER - COMM_OPS
    print(f"seam-grade floor at d={d}:  eps <= ({MULTIPLIER} - {COMM_OPS}) "
          f"p_loc = {share} p_loc")
    for name, p in P_LOC.items():
        eps = share * p
        print(f"  p_loc = {p:.0e}  ({name}):  eps_max = {eps:.1e}"
              f"  ->  F_raw >= {1-eps:.3f}")
    # asserted against the 5.3 prose
    assert share == 2
    assert abs(share * 1e-3 - 2e-3) < 1e-12
    assert abs(share * 3e-3 - 6e-3) < 1e-12

    print(f"\nmemory line, full T_round/T2 (stored half waits one round):")
    for k, br in enumerate(BRACKETS):
        tr = T.merge_window_us(d, k) * 1e-6
        line = " ".join(f"{tr/t2:.1e} ({nm.split(',')[0]})"
                        for nm, t2 in T2_S.items())
        print(f"  {br:14s} window {tr*1e3:5.1f} ms:  {line}")
    tr_base = T.merge_window_us(d, 1) * 1e-6
    assert 2e-4 < tr_base / 50.0 < 3e-4   # ~2.8e-4, the prose number
    floor_tight = share * 1e-3
    # The distilled lane's accumulation wait (thesis 4.3.4/7.3): the first-caught
    # half waits up to three windows before its batch fires, not one. The extra
    # two windows ride the RAW inputs, where the two selection checks catch
    # first-order errors, so only the survivor's single post-keep window reaches
    # the output; that is the one-round charge above. Checked here: the extra
    # input idle is small against the 6% raw error the distillation budget
    # already carries.
    extra_in = 2.0 * (T.merge_window_us(D, 1) * 1e-6) / 50.0
    assert extra_in < 0.01 * 0.06, "accumulation wait must be negligible against the raw error budget"
    print(f"  distilled-lane accumulation wait, 2 extra windows on the raw inputs: {extra_in:.1e},")
    print(f"  {extra_in/0.06*100:.1f}% of the 6% raw error the x3.4 demand already carries; first-order")
    print(f"  errors are caught by the selection checks, so the output keeps the one-round charge.")
    print(f"\nwait as a fraction of the tight floor {floor_tight:.0e}:")
    for k, br in enumerate(BRACKETS):
        tr = T.merge_window_us(d, k) * 1e-6
        print(f"  {br:14s} {tr/50.0:.1e}  = {100*(tr/50.0)/floor_tight:4.1f}% of floor")
    print("still charged in full and non-binding, but not negligible at the")
    print("conservative round; the prose states the fractions honestly (G2).")

    # The idle-dephasing charge (Section 4.2). Beyond the stored half's full-round wait, the
    # schedule leaves a cross-row or comm ancilla elevated in the transport lane for up to
    # worst_elevated_idle_us before its gate fires. The 137Ba+ clock qubit is first-order
    # field-insensitive, so that wait dephases at the same memory T2 whether the ion sits in a
    # memory well or elevated in the gradient. Charged in full at the tight 50 s T2, both idle
    # terms stay under the seam floor, so the raw visibility floor is set by the seam grade.
    print(f"\nidle dephasing, worst ion, charged in full at the clock T2 = 50 s:")
    for k, br in enumerate(BRACKETS):
        t_stored = T.merge_window_us(d, k) * 1e-6
        t_elev = T.worst_elevated_idle_us(d, k) * 1e-6
        e_stored, e_elev = t_stored / 50.0, t_elev / 50.0
        print(f"  {br:14s} stored half {t_stored*1e3:5.1f} ms -> {e_stored:.1e}, "
              f"elevated ancilla {t_elev*1e3:4.2f} ms -> {e_elev:.1e}  "
              f"(sum {e_stored+e_elev:.1e})")
    e_idle_base = (T.merge_window_us(d, 1) * 1e-6
                   + T.worst_elevated_idle_us(d, 1) * 1e-6) / 50.0
    assert e_idle_base < share * 1e-3, "idle dephasing must stay under the seam floor"
    print(f"sum at baseline {e_idle_base:.1e} stays under the {share*1e-3:.0e} seam floor, so the")
    print("raw visibility floor is set by the seam grade, not by idle dephasing.")

    # required visibility: F = (1+V)/2 when distinguishability is the only
    # imperfection, so eps = (1-V)/2 and the floor inverts to V >= 1-2 eps.
    print("\nrequired two-photon visibility, V >= 1 - 2 eps_max:")
    for name, p in P_LOC.items():
        v = 1.0 - 2.0 * share * p
        print(f"  p_loc = {p:.0e}  ({name}):  V >= {100*v:.1f}%")
    assert abs((1 - 2*share*1e-3) - 0.996) < 1e-12
    assert abs((1 - 2*share*3e-3) - 0.988) < 1e-12

    # The distilled operating floor, the point the machine runs at. Raw pairs feeding one
    # distillation round need only the simple-protocol tolerance, about ten times the local
    # gate error (MULTIPLIER), so the visibility requirement relaxes from the raw floor to
    # V >= 1 - 2 (MULTIPLIER p_loc). At the tight anchor this is the 98% of Section 5.4.
    print("\ndistilled operating floor, V >= 1 - 2 (MULTIPLIER p_loc):")
    for name, p in P_LOC.items():
        vd = 1.0 - 2.0 * MULTIPLIER * p
        print(f"  p_loc = {p:.0e}  ({name}):  V >= {100*vd:.1f}%")
    assert abs((1 - 2*MULTIPLIER*1e-3) - 0.98) < 1e-12   # the 98% distilled operating floor

    # The mismatch bound (Section 7.2). The circuit-level injection runs on the
    # abstract code, so the physical schedule's transport and idle steps are
    # unpriced there. Fold them into an effective local error
    # p_eff = p_loc + Delta. The distilled floor is not a cliff; it slides as
    # V >= 1 - 2 MULTIPLIER p_eff. The priced Delta is the worst-ion idle plus
    # the expected stalled-window charge (0.06 windows/round, make_requirement.py).
    p_anchor = 1e-3
    delta_priced = e_idle_base * 1.06
    v_slid = 1 - 2 * MULTIPLIER * (p_anchor + delta_priced)
    print("\nmismatch bound, p_eff = p_loc + Delta (schedule steps the abstract code omits):")
    print(f"  distilled floor slides as V >= 1 - {2*MULTIPLIER} p_eff")
    print(f"  priced Delta at baseline: idle {e_idle_base:.1e} + stall {0.06*e_idle_base:.1e}")
    print(f"  at the tight anchor {p_anchor:.0e} the floor slides 98.0% -> {100*v_slid:.1f}%")
    print(f"  at the family best gate 3e-4 the anchor absorbs Delta <= {p_anchor-3e-4:.1e}"
          f"  (x{(p_anchor-3e-4)/e_idle_base:.1f} the priced idle)")
    assert abs((1 - 2*MULTIPLIER*p_anchor) - 0.98) < 1e-12
    assert 0.972 < v_slid < 0.975

    # chain fidelity charges at the demonstrated link (O'Reilly et al. 2024,
    # PRL 133, 090802): polarization mixing dominated the pair error;
    # temporal mismatch and dark counts were bundled at 0.4%. Cited in 5.1.
    mixing, temporal_dark = 0.029, 0.004
    print(f"\nchain charges at the demonstrated link: "
          f"mixing {mixing:.1%}, temporal+dark {temporal_dark:.1%}")
    print(f"together {mixing+temporal_dark:.1%} vs the {share*1e-3:.1%} floor:")
    print("the demonstrated chain alone overspends the tight floor, so the")
    print("interface must beat demonstrated polarization handling, or the")
    print("link must distill (Section 5.4).")

"""
Chip-geometry frame generator + visualizer, driven by the scheduler.
====================================================================
Every page draws the full Chapter-4 cell geometry at its own positions:
the 2d+1 memory homes, the gate strip (junction columns and gate wells
alternating, then the three hold wells and the swap well), d SPAM sites,
the optical wall, and the cavity with its Yb coolant, one cell per code row.

Three families, at every odd distance:
  * round         -> one local syndrome-extraction round.
  * merge_full    -> two rounds of a remote lattice-surgery merge (raw lane).
  * merge_distill -> the same merge on the distilled lane (Section 4.3.4):
                     two halves ferried to the hold wells during the round,
                     the spent survivor recycled at the read to catch the
                     third, double selection at the boundary.

The scheduler's frame_errors legality check runs on every frame, at build time
and precomputed into each page's live line, so the page reflects the scheduler
and adds no rule of its own: well occupancy within capacity, no rest on a
junction column, and no reordering along a row without a shared-well rotation.
Run:  python3 qec_visualizer.py            # d = 3 and 5
      python3 qec_visualizer.py 7          # one distance
      python3 qec_visualizer.py all        # the thesis sweep, 3..27
"""
from __future__ import annotations
import json, math, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from qec_scheduler import build_stabilizers, place, stab_cell, num, seam_schedule, round_ops, parallel_steps
import qec_scheduler as S

MERGE_ROUNDS = 2
NAMES3 = {frozenset({1,2,4,5}):"A1", frozenset({2,3,5,6}):"A2",
          frozenset({4,5,7,8}):"A3", frozenset({5,6,8,9}):"A4",
          frozenset({2,3}):"B1", frozenset({7,8}):"B2",
          frozenset({1,4}):"B3", frozenset({6,9}):"B4"}

def make_labels(stabs, d):
    if d == 3:
        return {s: NAMES3[frozenset(num(rc, d) for rc in s.data)] for s in stabs}
    out, a, b = {}, 1, 1
    for s in stabs:
        if s.weight == 4: out[s] = "A%d" % a; a += 1
        else:             out[s] = "B%d" % b; b += 1
    return out

P = 46.0
IW = 9.0

class G:
    """The chip geometry at one distance, in pixels."""
    def __init__(self, d):
        self.d = d
        self.MX = [80 + k*P for k in range(2*d + 1)]
        x = self.MX[-1] + 74
        self.JCOL, self.WELL = {}, {}
        for c in range(d):
            self.JCOL[c] = x; x += P
            self.WELL[c] = x; x += P
        self.HOLDX = [x, x + P, x + 2*P]; x += 3*P
        self.SWX = x + 6; x += P + 12
        npark = max(1, (d - 1)//2)
        self.PARKX = [x + k*P for k in range(npark)]; x += npark*P + 20
        self.SPX = [x + k*P for k in range(d)]; x += d*P + 30
        self.WALLX = x; x += 44
        self.CAVX = x + 12; self.YBX = self.CAVX + 32
        self.CY = {r: 110 + 190*r for r in range(d)}
        self.SLOTS = (self.MX + list(self.WELL.values()) + self.HOLDX + [self.SWX]
                      + self.PARKX + self.SPX + [self.CAVX, self.YBX])
        self.CAP = {**{x2: 2 for x2 in self.MX},
                    **{x2: 3 for x2 in self.WELL.values()},
                    **{x2: 2 for x2 in self.HOLDX}, self.SWX: 2,
                    **{x2: 2 for x2 in self.PARKX},
                    **{x2: 2 for x2 in self.SPX}, self.CAVX: 2, self.YBX: 1}
    def gapy(self, a, b): return (self.CY[a] + self.CY[b]) / 2
    def rule_geom(self):
        """The pixel descriptor the scheduler's frame_errors reads, so the
        physical rules live in the scheduler and this class only supplies
        coordinates."""
        return {"d": self.d, "JCOL": list(self.JCOL.values()),
                "CY": {r: self.CY[r] for r in range(self.d)}, "CAP": dict(self.CAP),
                "WELL": list(self.WELL.values()), "SPX": list(self.SPX),
                "SWX": self.SWX, "WALLX": self.WALLX}


CHECK_ONLY = False

def build(d, merge, rounds, distill=False):
    STABS = build_stabilizers(d)
    LABEL = make_labels(STABS, d)
    SEAM = seam_schedule(d)
    CHAINS = place(d)
    ID_DATA = {(r, c): "d%d" % (r*d + c + 1) for r in range(d) for c in range(d)}
    def A(s): return LABEL[s]
    def Dt(rc): return ID_DATA[rc]
    g = G(d)
    CY = g.CY

    ions, home = {}, {}
    for ci, chain in enumerate(CHAINS):
        for j, (kind, item) in enumerate(chain):
            i = ID_DATA[item] if kind == "data" else LABEL[item]
            ions[i] = "data" if kind == "data" else ("X" if item.kind == "X" else "Z")
            if kind == "data":
                home[i] = (g.MX[1 + 2*item[1]], ci)
            else:
                lft = chain[j - 1] if j > 0 else None
                rgt = chain[j + 1] if j < len(chain) - 1 else None
                if lft and lft[0] == "data" and rgt and rgt[0] == "data":
                    home[i] = (g.MX[2 + 2*lft[1][1]], ci)
                elif lft is None:
                    home[i] = (g.MX[0], ci)
                else:
                    home[i] = (g.MX[2*d], ci)
    for r in range(d):
        ions["C%d" % r] = "comm" if r in SEAM else "spare"
        home["C%d" % r] = (g.CAVX, r)
        if merge and r in SEAM:
            ions["C%dB" % r] = "comm"; home["C%dB" % r] = (g.SWX, r)
            if distill:
                for k, hx in enumerate(g.HOLDX):
                    hid = "h%d" % (3*r + k + 1)
                    ions[hid] = "held"; home[hid] = (hx, r)
        ions["Y%d" % r] = "yb"; home["Y%d" % r] = (g.YBX, r)

    slot = {i: home[i] for i in home}
    FR = []
    ROT = []
    TRANSIT = set()

    def groups(row):
        gg = {}
        for i, s in slot.items():
            if len(s) == 2 and not isinstance(s[0], str) and s[1] == row:
                gg.setdefault(s[0], []).append(i)
        return gg

    def render():
        posr = {}
        gg = {r: groups(r) for r in range(d)}
        for i, s in slot.items():
            if len(s) == 3:
                posr[i] = [s[1], s[2]]
            else:
                x, row = s
                mem = sorted(gg[row].get(x, [i]))
                n = len(mem); k = mem.index(i)
                off = 0 if n == 1 else (k - (n - 1)/2) * 2*IW
                posr[i] = [round(x + off, 1), CY[row]]
        return posr

    ERRS = []
    RULEGEOM = g.rule_geom()
    _last_frame = [None]

    def _inc_check(op=None, pairs=None, tv=None):
        # Reflect the scheduler's rules on the current frame (streaming, so the
        # big sweep keeps only the previous frame in memory).
        fs = {i: (list(slot[i]) if len(slot[i]) == 3 else [slot[i][0], slot[i][1]]) for i in slot}
        f = {"slots": fs}
        if op: f["op"] = op
        if pairs: f["pairs"] = [sorted(pp) for pp in pairs]
        tvs = sorted(TRANSIT | ({tv} if tv else set()))
        if tvs: f["tv"] = tvs
        fi = len(FR)
        for m in S.per_frame_errors(f, RULEGEOM, ions):
            ERRS.append((fi, m))
        if _last_frame[0] is not None:
            for m in S.pair_errors(_last_frame[0], f, RULEGEOM):
                ERRS.append((fi, m))
        _last_frame[0] = f

    def snap(cap, hi=None, junc=None, badge="", op=None, pairs=None, tv=None):
        if CHECK_ONLY:
            _inc_check(op, pairs, tv)
            FR.append(1)
            return
        f = {"slots": {i: (list(slot[i]) if len(slot[i]) == 3 else [slot[i][0], slot[i][1]]) for i in slot},
             "pos": render(), "cap": cap, "hi": hi or [], "junc": junc or [], "badge": badge}
        if op: f["op"] = op
        if pairs: f["pairs"] = [sorted(p) for p in pairs]
        tvs = sorted(TRANSIT | ({tv} if tv else set()))
        if tvs: f["tv"] = tvs
        FR.append(f)

    def occ(x, row):
        return [i for i, s in slot.items() if len(s) == 2 and s[0] == x and s[1] == row]

    def hop(i, tx, row, cap, hi=None, badge=""):
        while len(occ(tx, row)) >= g.CAP.get(tx, 99) and i not in occ(tx, row):
            e = [x2 for x2 in occ(tx, row) if x2 != i][0]
            region = (list(g.WELL.values()) if tx in g.WELL.values() else
                      g.MX if tx in g.MX else
                      g.HOLDX + [g.SWX] if (tx in g.HOLDX or tx == g.SWX) else
                      g.SPX if tx in g.SPX else g.SLOTS)
            spots = sorted([w for w in region if w != tx and len(occ(w, row)) < g.CAP.get(w, 0)],
                           key=lambda w: abs(w - tx))
            hop(e, spots[0], row, cap + f" {e} steps aside within its zone to make room.", hi=[e])
        TRANSIT.add(i)
        x0 = slot[i][0] if len(slot[i]) == 2 else slot[i][1]
        step = 1 if tx > x0 else -1
        path = sorted([x for x in g.SLOTS if (x - x0)*step > 0.5 and (tx - x)*step > 0.5],
                      key=lambda x: (x - x0)*step)
        for bx in path:
            if occ(bx, row):
                slot[i] = (bx, row)
                snap(cap + f" {i} merges into the occupied well on its path and pairwise rotations carry it through.",
                     hi=hi or [i], badge=badge or "rotate-through", tv=i)
                if bx == g.SWX: ROT.append(len(FR) - 1)
        slot[i] = (tx, row)
        TRANSIT.discard(i)
        snap(cap, hi=hi or [i], badge=badge)

    def ensure_room(tx, row, cap):
        while len(occ(tx, row)) >= g.CAP.get(tx, 99):
            e = occ(tx, row)[0]
            region = (list(g.WELL.values()) if tx in g.WELL.values() else
                      g.MX if tx in g.MX else
                      g.HOLDX + [g.SWX] if (tx in g.HOLDX or tx == g.SWX) else
                      g.SPX if tx in g.SPX else g.SLOTS)
            spots = sorted([w for w in region if w != tx and len(occ(w, row)) < g.CAP.get(w, 0)],
                           key=lambda w: abs(w - tx))
            hop(e, spots[0], row, cap + f" {e} steps aside within its zone to make room.", hi=[e])

    def multi(assign, cap, hi=None, badge=""):
        for i, (x, row) in assign.items(): slot[i] = (x, row)
        snap(cap, hi=hi or list(assign), badge=badge)

    def hop_to_channel(i, jx, row, cap, hi=None):
        """Travel along the row to the junction column and pause in the channel
        there, rotating through every occupied well on the way."""
        x0 = slot[i][0] if len(slot[i]) == 2 else slot[i][1]
        step = 1 if jx > x0 else -1
        path = sorted([x for x in g.SLOTS if (x - x0)*step > 0.5 and (jx - x)*step > 0.5],
                      key=lambda x: (x - x0)*step)
        TRANSIT.add(i)
        for bx in path:
            if occ(bx, row):
                slot[i] = (bx, row)
                snap(cap + f" {i} merges into the occupied well on its path and pairwise rotations carry it through.",
                     hi=hi or [i], badge="rotate-through", tv=i)
                if bx == g.SWX: ROT.append(len(FR) - 1)
        slot[i] = ("J", jx, CY[row])
        TRANSIT.discard(i)
        snap(cap, hi=hi or [i], badge="at the junction column")

    def cross_rows(i, c, src, dst, cap, hi=None):
        """A full junction transit: along the source row to junction column c,
        vertically through the junction, landing in the destination row's
        channel at the same column."""
        hop_to_channel(i, g.JCOL[c], src, cap + f" {i} moves along its row to junction column {c}.", hi=hi)
        slot[i] = ("J", g.JCOL[c], g.gapy(min(src, dst), max(src, dst)))
        snap(cap + f" {i} lifts through the junction between the rows.",
             hi=hi or [i], junc=[[c, min(src, dst)]], badge="in transit")
        slot[i] = ("J", g.JCOL[c], CY[dst])
        snap(cap + f" {i} arrives in the destination row's channel at the junction column.",
             hi=hi or [i], badge="in transit")

    BANDW = {}
    def band_in():
        BANDW.clear()
        assign = {}
        for row in range(d):
            chain = sorted([i for i in slot if len(slot[i]) == 2 and slot[i][1] == row
                            and slot[i][0] in g.MX and ions[i] in ("data", "X", "Z")],
                           key=lambda i: slot[i][0])
            gs = [chain[k:k+2] for k in range(0, len(chain), 2)]
            for w, grp in enumerate(gs):
                for i in grp:
                    assign[i] = (g.WELL[min(w, d - 1)], row); BANDW[i] = g.WELL[min(w, d - 1)]
        multi(assign, "Gate feed: each row's band files into the gate strip in order, two ions to a well. One entry, as priced.",
              badge="band enters strip")

    def band_out():
        movers = sorted([i for i in BANDW if slot[i] != home[i]],
                        key=lambda i: home[i][0])
        for i in movers:
            hop(i, home[i][0], home[i][1], f"Band return: {i} files back to its memory home.", hi=[i])
        snap("The band is back in its memory homes.", hi=movers, badge="band returns")

    def to_well(i, wx, row, cap, hi=None):
        if slot[i] != (wx, row):
            hop(i, wx, row, cap, hi=hi)

    def isolate(wx, row, keep, avoid, cap):
        """Evict every resident of the well except `keep`, each to the nearest
        well with room, preferring wells not gating this step, so the gate
        fires on the pair alone."""
        extras = [i for i in occ(wx, row) if i not in keep]
        for e in extras:
            spots = sorted([w for w in g.WELL.values()
                            if w != wx and len(occ(w, row)) < g.CAP[w]],
                           key=lambda w: (w in avoid, abs(w - wx)))
            hop(e, spots[0], row, cap + f" {e} steps aside to a well with room so the pair is alone.", hi=[e])

    def gate_pair(a, dd, row, avoid, cap):
        wx = slot[a][0] if slot[a][0] in g.WELL.values() else slot[dd][0]
        isolate(wx, row, {a, dd}, avoid, cap)
        if slot[dd] != (wx, row):
            hop(dd, wx, row, cap + f" {dd} steps into the gate well.", hi=[dd])
        if slot[a] != (wx, row):
            hop(a, wx, row, cap + f" {a} steps into the gate well.", hi=[a])
        return wx

    net = {l: "C%d" % l for l in SEAM}
    car = {l: "C%dB" % l for l in SEAM}
    surv, PEND, FILEOUT, XPEND = {}, {}, {}, {}
    FERRYQ, NEED3, CHK = {}, {}, {}
    def C(l): return surv[l] if (distill and l in surv) else car[l]

    def free_hold_bodies(l):
        return sorted([i for i in ions if ions[i] == "held" and slot[i][1] == l and slot[i][0] in g.HOLDX],
                      key=lambda i: slot[i][0])

    def ferry_half(l, k, hh, tag):
        snap(f"{tag} lane {l}: the networker heralds raw pair {k} at the cavity.", hi=[net[l]], badge="herald")
        hop(net[l], g.SWX, l, "The networker carries the fresh half to the swap well.")
        snap(f"Crystal rotation at the swap well: the carrier {car[l]} takes the half.", hi=[net[l], car[l]])
        hop(net[l], g.CAVX, l, "The networker returns to the cavity and keeps attempting.")
        hop(car[l], slot[hh][0], l, f"The carrier parks half {k} beside {hh} in a gate-end hold well, a one-hop ferry.")
        hop(car[l], g.SWX, l, "The carrier returns to the swap well.")

    def distill_cnots(l, order, tag):
        sv, c1, c2 = order
        hop(sv, g.WELL[1], l, f"{tag} lane {l}: {sv}, the first catch, takes the survivor's middle gate well.")
        hop(c1, g.WELL[0], l, f"{c1} takes the left gate well.")
        hop(c2, g.WELL[2], l, f"{c2} takes the right gate well.")
        slot[c1] = (g.WELL[1], l)
        snap(f"Bilateral CNOT one: {c1} and the survivor, an isolated pair.", hi=[c1, sv],
             op="gate", pairs=[(c1, sv)])
        slot[c1] = (g.WELL[0], l)
        snap(f"{c1} splits back out.", hi=[c1])
        slot[c2] = (g.WELL[1], l)
        snap(f"Bilateral CNOT two: {c2} and the survivor, an isolated pair.", hi=[c2, sv],
             op="gate", pairs=[(c2, sv)])
        slot[c2] = (g.WELL[2], l)
        snap(f"{c2} splits back out and takes its basis rotation.", hi=[c2])
        surv[l] = sv

    def distill_reads(l, c1, c2, sv):
        hop(c1, g.SPX[0], l, f"{c1} hops to a SPAM site and is read in Z, the bit-flip check.")
        snap(f"Read of {c1} in Z.", hi=[c1], badge="read Z", op="read", pairs=[(c1,)])
        hop(c2, g.SPX[1], l, f"{c2} follows and is read in X, the phase-flip check.")
        snap(f"Read of {c2} in X.", hi=[c2], badge="read X", op="read", pairs=[(c2,)])
        snap(f"Both modules' checks agree: keep. {sv} holds the purified pair.", hi=[sv], badge="keep")
        hop(c1, g.HOLDX[1], l, f"{c1} resets and returns to the hold pool.")
        hop(c2, g.HOLDX[2], l, f"{c2} resets and returns to the hold pool.")

    _ops = round_ops(d, merge, rounds)
    for op in _ops:
        v = op[0]
        if v == "prep":
            snap("Rest. Memory homes hold the code row"
                 + (", the carrier waits at each active swap well, the networker and Yb sit at each cavity." if merge
                    else ". The comm ions rest at their cavities; a local round asks nothing of them.")
                 + (" The three hold wells beside each swap well serve the distilling lanes." if distill else ""),
                 badge="chapter-4 geometry")
        elif v == "round":
            if distill and not surv:
                snap("Warm-up: with the blocked lanes cleared, each active lane distills its first batch before round 1.", badge="warm-up")
                for l in sorted(SEAM):
                    pool = free_hold_bodies(l)
                    for k, hh in enumerate(pool):
                        ferry_half(l, k, hh, "Warm-up")
                    distill_cnots(l, pool, "Warm-up")
                    distill_reads(l, pool[1], pool[2], pool[0])
                    hop(pool[0], g.SWX, l, f"{pool[0]} joins the carrier at the swap well, purified pair ready for the seam.")
            snap(f"Round {op[1] + 1} of {op[2]}.", badge=f"round {op[1] + 1}")
        elif v == "park":
            for s, cl, well in op[1]:
                for row2 in range(cl, d - 1):
                    cross_rows(A(s), d - 1, row2, row2 + 1, f"Park:", hi=[A(s)])
                hop(A(s), g.PARKX[well], d - 1, f"{A(s)} moves along the bottom row into its spare park well.", hi=[A(s)])
        elif v == "unpark":
            for s, cl, well in reversed(op[1]):
                for row2 in range(d - 1, cl, -1):
                    cross_rows(A(s), d - 1, row2, row2 - 1, f"Unpark:", hi=[A(s)])
                hop(A(s), home[A(s)][0], cl, f"{A(s)} walks along its home row back to its memory well.", hi=[A(s)])

        elif v in ("inrow","swap","xlift","xgate","xlower","xdrop","comm_out","comm_lift",
                   "comm_gate","comm_lower","comm_back","comm_arrive"):
            if v == "inrow" and not BANDW:
                band_in()
            if v == "inrow":
                L = op[1] + 1
                avoid = {slot[A(s)][0] for s, _ in op[2] if len(slot[A(s)]) == 2}
                for s, rc in op[2]:
                    gate_pair(A(s), Dt(rc), slot[A(s)][1], avoid, f"Step {L}:")
                    snap(f"Step {L}: the in-row gate fires, an isolated two-ion crystal in a gate well.",
                         hi=[A(s), Dt(rc)], op="gate", pairs=[(A(s), Dt(rc))])
            elif v == "swap":
                L = op[1] + 1
                for s, rc in op[2]:
                    a, dd = A(s), Dt(rc)
                    src, dst = slot[a], slot[dd]
                    hop(dd, src[0], src[1], f"Step {L}: {dd} works its way over, rotating past {a}.", hi=[dd])
                    hop(a, dst[0], dst[1], f"Step {L}: {a} takes {dd}'s freed well; the two have exchanged.", hi=[a])
            elif v == "xlift":
                L = op[1] + 1
                XPEND.clear()
                for s, c, ac, tr in op[2]:
                    hop_to_channel(A(s), g.JCOL[c], ac, f"Step {L}:", hi=[A(s)])
                    slot[A(s)] = ("J", g.JCOL[c], g.gapy(min(ac, tr), max(ac, tr)))
                    XPEND[A(s)] = (c, ac, tr)
                snap(f"Step {L}: the cross-row ancillas hang in their junctions at c + 1/4, between the rows.",
                     hi=[A(s) for s, *_ in op[2]], junc=[[c, min(ac, tr)] for s, c, ac, tr in op[2]], badge="in transit")
            elif v == "xgate":
                L = op[1] + 1
                for s, rc in op[2]:
                    c, ac, tr = XPEND[A(s)]
                    isolate(g.WELL[c], rc[0], {Dt(rc), A(s)}, set(), f"Step {L}:")
                    if slot[Dt(rc)][0] != g.WELL[c]:
                        hop(Dt(rc), g.WELL[c], rc[0], f"Step {L}: {Dt(rc)} steps into the gate well beside the junction.", hi=[Dt(rc)])
                    slot[A(s)] = ("J", g.JCOL[c], CY[rc[0]])
                    snap(f"Step {L}: {A(s)} drops into the prepared row's channel at its junction column.",
                         hi=[A(s)], badge="in transit")
                    hop(A(s), g.WELL[c], rc[0], f"Step {L}: {A(s)} joins {Dt(rc)}; the pair is alone in the well.", hi=[A(s)])
                    snap(f"Step {L}: the cross-row gate fires on the isolated pair.",
                         hi=[A(s)], op="gate", pairs=[(A(s), Dt(rc))])
            elif v == "xlower":
                L = op[1] + 1
                for s, c, ac, tr in op[2]:
                    hop_to_channel(A(s), g.JCOL[c], tr, f"Step {L}:", hi=[A(s)])
                    slot[A(s)] = ("J", g.JCOL[c], g.gapy(min(ac, tr), max(ac, tr)))
                snap(f"Step {L}: they lift back into their junctions.",
                     hi=[A(s) for s, *_ in op[2]], junc=[[c, min(ac, tr)] for s, c, ac, tr in op[2]], badge="in transit")
            elif v == "xdrop":
                L = op[1] + 1
                for s, c, ac, swapped in op[2]:
                    tgt = BANDW.get(A(s), g.WELL[d - 1])
                    ensure_room(tgt, ac, f"Step {L}:")
                    slot[A(s)] = ("J", g.JCOL[c], CY[ac])
                    snap(f"Step {L}: {A(s)} drops back into its home row's channel.", hi=[A(s)], badge="in transit")
                    hop(A(s), tgt, ac, f"Step {L}: {A(s)} returns along the row to its own strip well.", hi=[A(s)])
            elif v == "comm_out":
                L = op[1] + 1
                snap(f"Step {L}: the seam ions ready at their swap wells while the seam gate wells clear.",
                     hi=[C(l) for l in op[2]], badge="seam excursion")
            elif v == "comm_lift":
                L = op[1] + 1
                for l in op[2]:
                    hop_to_channel(C(l), g.JCOL[d - 1], l, f"Step {L}:", hi=[C(l)])
                    slot[C(l)] = ("J", g.JCOL[d - 1], g.gapy(l, l + 1))
                snap(f"Step {L}: the seam ions hang in the boundary junctions, between the rows.",
                     hi=[C(l) for l in op[2]], junc=[[d - 1, l] for l in op[2]], badge="in transit")
            elif v == "comm_gate":
                L = op[1] + 1
                for l, rc in op[2]:
                    isolate(g.WELL[d - 1], rc[0], {Dt(rc), C(l)}, set(), f"Step {L}:")
                    if slot[Dt(rc)][0] != g.WELL[d - 1]:
                        hop(Dt(rc), g.WELL[d - 1], rc[0], f"Step {L}: the boundary data steps into the seam gate well.", hi=[Dt(rc)])
                    if len(slot[C(l)]) == 3 and rc[0] != l:
                        slot[C(l)] = ("J", g.JCOL[d - 1], CY[rc[0]])
                        snap(f"Step {L}: {C(l)} drops into the prepared row's channel at the boundary junction.",
                             hi=[C(l)], badge="in transit")
                    hop(C(l), g.WELL[d - 1], rc[0], f"Step {L}: {C(l)} joins the boundary data, the pair alone in the well.", hi=[C(l)])
                    snap(f"Step {L}: the seam gate fires on the isolated pair.",
                         hi=[C(l)], op="gate", pairs=[(C(l), Dt(rc))])
                if distill:
                    for l, rc in op[2]:
                        if FERRYQ.get(l):
                            FERRYQ[l] = False
                            pool = free_hold_bodies(l)[:2]
                            for k, hh in enumerate(pool):
                                ferry_half(l, k, hh, "With the seam ion out, next batch,")
                            PEND[l] = pool
            elif v == "comm_lower":
                L = op[1] + 1
                for l in op[2]:
                    hop_to_channel(C(l), g.JCOL[d - 1], l + 1, f"Step {L}:", hi=[C(l)])
                    slot[C(l)] = ("J", g.JCOL[d - 1], g.gapy(l, l + 1))
                snap(f"Step {L}: the seam ions lift back into the junctions.",
                     hi=[C(l) for l in op[2]], junc=[[d - 1, l] for l in op[2]], badge="in transit")
            elif v == "comm_back":
                L = op[1] + 1
                for l in op[2]:
                    if len(slot[C(l)]) == 3:
                        slot[C(l)] = ("J", g.JCOL[d - 1], CY[l])
                        snap(f"Step {L}: the seam ion drops back into its home row's channel.", hi=[C(l)], badge="in transit")
                    hop(C(l), g.SWX, l, f"Step {L}: it returns along the row to the gate-end swap well.", badge="seam excursion")
            elif v == "comm_arrive":
                for l in op[2]:
                    if slot[C(l)] != (g.SWX, l):
                        hop(C(l), g.SWX, l, "The seam ion settles at its swap well.")
        elif v == "readout":
            if BANDW:
                band_out(); BANDW.clear()
                if distill and NEED3:
                    for l in sorted(NEED3):
                        old = NEED3.pop(l)
                        ferry_half(l, 2, old, "Third catch, off the critical path:")
                        distill_cnots(l, PEND.get(l, [])[:2] + [old], "Boundary:")
                        CHK[l] = (PEND[l][0], PEND[l][1], old)
            for s, rc in op[1]:
                a, dd = A(s), Dt(rc)
                sa, sd = slot[a], slot[dd]
                hop(a, sd[0], sd[1], f"File-out: {a} merges into {dd}'s well.", hi=[a])
                hop(dd, sa[0], sa[1], f"The crystal rotates and splits: {a} has passed {dd}.", hi=[dd])
        elif v == "to_spam":
            order = sorted([A(s) for s in op[1]], key=lambda i: -slot[i][0])
            FILEOUT.update({a: slot[a] for a in order})
            for a in order:
                row = slot[a][1]
                free = [x for x in g.SPX if not occ(x, row)]
                hop(a, free[0], row, f"{a} shuttles the corridor, across the strip, to a SPAM site.", badge="convoy")
            snap("All ancillas sit at SPAM detection sites, behind the wall from the memory data.", hi=order, badge="in SPAM")
        elif v == "syndromes":
            snap("493 nm readout fires at the SPAM sites"
                 + (", the seam bits read at the swap wells in the same beat." if op[2] else "."),
                 hi=[A(s) for s in op[1]], badge="syndromes")
        elif v == "measure":
            if distill:
                snap("The survivors are read at the swap wells. This round's purified pairs reach module B.",
                     hi=[C(l) for l in op[1]], badge="purified pairs to B", op="read",
                     pairs=[(C(l),) for l in op[1]])
                for l in op[1]:
                    old = surv.pop(l, None)
                    if old is None:
                        continue
                    free = [x for x in g.HOLDX if not occ(x, l)]
                    hop(old, free[0], l, f"{old} resets and rejoins the hold pool, an empty body again.")
                    NEED3[l] = old
            else:
                snap("The carriers are read at their swap wells. This round's Bell pairs reach module B.",
                     hi=[C(l) for l in op[1]], badge="pairs to B", op="read",
                     pairs=[(C(l),) for l in op[1]])
        elif v == "from_spam":
            order = sorted([A(s) for s in op[1]], key=lambda a: -FILEOUT[a][0])
            for a in order:
                tx, row = FILEOUT[a]
                hop(a, tx, row, f"{a} shuttles back in from the SPAM zone to where it left the row.", badge="convoy")
            if distill and CHK:
                for l in sorted(CHK):
                    sv = surv[l]
                    trio = CHK.pop(l)
                    c1, c2 = trio[1], trio[2]
                    distill_reads(l, c1, c2, sv)
        elif v == "reset":
            for s, rc in op[1]:
                a, dd = A(s), Dt(rc)
                sa, sd = slot[a], slot[dd]
                hop(a, sd[0], sd[1], f"Reset file-in: {a} merges back into {dd}'s well.", hi=[a])
                hop(dd, sa[0], sa[1], f"Rotation undone: {a} is back on the memory side of {dd}.", hi=[dd])
        elif v == "reset_done":
            if distill:
                for l in sorted(SEAM):
                    if l in surv and slot[surv[l]] != (g.SWX, l):
                        hop(surv[l], g.SWX, l, f"{surv[l]} joins the carrier at the swap well, ready for the next seam round.")
            strays = [i for i in ions if ions[i] in ("data", "X", "Z")
                      and i in home and slot[i] != home[i]
                      and not (merge and any(slot[i] == (px, d - 1) for px in g.PARKX))]
            for i in sorted(strays, key=lambda i: home[i][0]):
                hop(i, home[i][0], home[i][1], f"{i} walks back to its memory home.", hi=[i])
            snap("Everyone is back in the memory home it started from. The round ends at rest.", badge="at rest")
        elif v == "herald":
            if distill:
                for l in sorted(op[1]):
                    FERRYQ[l] = True
                snap("The networkers attempt through the round; the catches land once the seam ion "
                     "leaves the swap well, so the ferry rotations stay pairwise.",
                     hi=[net[l] for l in op[1]], badge="herald")
            else:
                snap("The networkers herald fresh Bell pairs at their cavities; Yb keeps them cold.",
                     hi=[net[l] for l in op[1]], badge="herald")
        elif v == "comm_swap":
            if distill:
                snap("Ping-pong roles hold in the distilled lane: the networker keeps the cavity, "
                     "the carrier keeps the ferry, and the survivor pipeline replaces the raw handoff.",
                     badge="distilled lane")
                continue
            for l in op[1]:
                hop(net[l], g.SWX, l, "Ping-pong: the heralded networker crosses the wall and the SPAM row to the swap well.", badge="handoff")
                snap("The two comm ions share the swap well and the crystal rotates: roles exchange.", hi=[net[l], car[l]])
                hop(car[l], g.CAVX, l, "The emptied ion returns to the cavity as the next networker.")
                net[l], car[l] = car[l], net[l]
        else:
            raise ValueError(f"no renderer for scheduler op {v!r}")
    build.last_rot = len(ROT)
    build.last_errs = ERRS
    return FR, ions, g


TEMPLATE = r"""<!doctype html><html lang="en"><head><meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>__TITLE__</title><style>
body{background:#232320;color:#ddd;font:14px -apple-system,'Segoe UI',sans-serif;margin:0;padding:14px}
#hdr b{font-size:16px}#hdr{margin:4px 8px 10px}.sub{color:#9a9a90;font-size:13px}
#stage{position:relative;background:#232320;border:1px solid #3a3a36;border-radius:10px;overflow:auto;max-height:78vh}
#world{position:relative}
.zone{position:absolute;border-radius:8px;border:1px dashed #4c4c46}
.zlbl{position:absolute;color:#8f8f86;font-size:12px;white-space:nowrap}
.well{position:absolute;width:34px;height:34px;border:1.4px solid #56564e;border-radius:9px;transform:translate(-50%,-50%)}
.well.gate{border-color:#4a6f96}.well.swapw{border-color:#7d74c9;border-style:dashed}
.well.hold{border-color:#3f8f71;border-style:dashed}
.well.spam{border-color:#a2793c}.well.cav{border-color:#3f8f71}.well.park{border-color:#777}
.jcol{position:absolute;width:0;border-left:2px dotted #5c5c54;transform:translateX(-50%)}
.jlink{position:absolute;width:2px;background:#4c4c46;transform:translateX(-50%)}
.jdot{position:absolute;width:7px;height:7px;border-radius:50%;background:#8a8a80;transform:translate(-50%,-50%)}
.wall{position:absolute;width:7px;background:#6e6e66;border-radius:2px;transform:translateX(-50%)}
.ion{position:absolute;width:27px;height:27px;border-radius:8px;display:flex;align-items:center;justify-content:center;
 font-weight:600;font-size:10.5px;transform:translate(-50%,-50%);transition:left .22s ease,top .22s ease;z-index:5;color:#fff}
.data{background:#2e2e2a;border:1.6px solid #9a9a90;color:#e8e8e0}
.X{background:#3f7fd4;border:1.6px solid #6fa5ec}.Z{background:#c96f3b;border:1.6px solid #e59a6c}
.comm{background:#2f8f6b;border:1.6px solid #5cc39a}.held{background:#26463a;border:1.6px solid #5cc39a;color:#bfe8d4}
.spare{background:#3a3a36;border:1.6px solid #6e6e66;color:#aaa}
.yb{width:16px;height:16px;border-radius:50%;background:#d85a30;border:1.4px solid #f08a60;font-size:0}
.hi{box-shadow:0 0 0 3px rgba(240,220,120,.85),0 0 14px rgba(240,220,120,.5);z-index:8}
#cap{min-height:44px;margin:12px 8px 6px;font-size:15px;line-height:1.45}
#badge{display:inline-block;background:#3a3a2f;border:1px solid #6b6b4f;color:#e6d98a;border-radius:20px;padding:2px 12px;font-size:12px;margin-left:10px;vertical-align:2px}
#verify{margin:0 8px;font-size:13px}.ok{color:#7ec98f}.bad{color:#e8806a}
#bar{display:flex;gap:10px;align-items:center;margin:10px 8px}
button{background:#32322e;color:#ddd;border:1px solid #55554d;border-radius:8px;padding:6px 14px;font-size:14px;cursor:pointer}
button:hover{background:#3d3d38}input[type=range]{flex:1}
.legend{margin:8px;color:#9a9a90;font-size:12.5px}
.lg{display:inline-block;width:13px;height:13px;border-radius:4px;vertical-align:-2px;margin:0 4px 0 12px}
</style></head><body>
<div id="hdr"><b>__TITLE__</b><br><span class="sub">__SUB__
The scheduler's own legality check is shown live on every frame: well capacity, junction columns clear of resting ions, and no passing without a shared-well crystal rotation.</span></div>
<div id="stage"><div id="world"></div></div>
<div id="cap"></div><div id="verify"></div>
<div id="bar"><button id="b0">&#9198;</button><button id="bp">&#9664;</button><button id="pl">Play</button>
<button id="bn">&#9654;</button><button id="be">&#9197;</button>
<input type="range" id="sl" min="0" value="0"><span id="ctr"></span>
<select id="sp"><option value="900">slow</option><option value="450" selected>normal</option><option value="180">fast</option></select></div>
<div class="legend"><span class="lg" style="background:#2e2e2a;border:1.6px solid #9a9a90"></span>data
<span class="lg" style="background:#3f7fd4"></span>X check <span class="lg" style="background:#c96f3b"></span>Z check
<span class="lg" style="background:#2f8f6b"></span>comm Ba <span class="lg" style="background:#26463a;border:1.6px solid #5cc39a"></span>held half
<span class="lg" style="background:#d85a30;border-radius:50%"></span>Yb coolant
<span class="lg" style="border:1.4px dashed #7d74c9;background:none"></span>swap well
<span class="lg" style="border:1.4px dashed #3f8f71;background:none"></span>hold wells
<span class="lg" style="border-left:2px dotted #5c5c54;background:none;width:2px"></span>junction column &middot;
memory homes 2d+1 &middot; strip: junction, well alternating, then holds and swap &middot; SPAM d sites &middot; wall &middot; cavity</div>
<script>
const G=__GEOM__, FR=__FRAMES__, KIND=__IONS__, CAP=__CAP__, VERDICTS=__VERDICTS__;
const W=document.getElementById('world');
const CY=Object.keys(G.CY).map(k=>G.CY[k]);
W.style.width=(G.YBX+80)+'px';W.style.height=(CY[CY.length-1]+90)+'px';
function zone(x1,x2,y,lbl){const z=document.createElement('div');z.className='zone';
 z.style.left=x1+'px';z.style.width=(x2-x1)+'px';z.style.top=(y-34)+'px';z.style.height='68px';W.appendChild(z);
 if(lbl){const t=document.createElement('div');t.className='zlbl';t.textContent=lbl;t.style.left=x1+'px';t.style.top=(y-52)+'px';W.appendChild(t);}}
CY.forEach((y,r)=>{
 zone(G.MX[0]-26,G.MX[G.MX.length-1]+26,y, r==0?'Memory (2d+1 homes)':'');
 zone(G.JCOL[0]-22,G.WELL[G.WELL.length-1]+24,y, r==0?'Gate strip':'');
 if(G.HOLDX.length)zone(G.HOLDX[0]-24,G.HOLDX[2]+24,y, r==0?'holds':'');
 zone(G.SWX-24,G.SWX+24,y, r==0?'swap':'');
 zone(G.SPX[0]-26,G.SPX[G.SPX.length-1]+26,y, r==0?'SPAM (d sites)':'');
 zone(G.CAVX-26,G.YBX+22,y, r==0?'Optical I/F':'');
 for(const x of G.MX){const w=document.createElement('div');w.className='well';w.style.left=x+'px';w.style.top=y+'px';W.appendChild(w);}
 for(const x of G.WELL){const w=document.createElement('div');w.className='well gate';w.style.left=x+'px';w.style.top=y+'px';W.appendChild(w);}
 for(const x of G.JCOL){const j=document.createElement('div');j.className='jcol';j.style.left=x+'px';j.style.top=(y-26)+'px';j.style.height='52px';W.appendChild(j);}
 for(const x of G.HOLDX){const w=document.createElement('div');w.className='well hold';w.style.left=x+'px';w.style.top=y+'px';W.appendChild(w);}
 {const w=document.createElement('div');w.className='well swapw';w.style.left=G.SWX+'px';w.style.top=y+'px';W.appendChild(w);}
 if(r==CY.length-1){for(const x of G.PARKX){const w=document.createElement('div');w.className='well park';w.style.left=x+'px';w.style.top=y+'px';W.appendChild(w);}}
 for(const x of G.SPX){const w=document.createElement('div');w.className='well spam';w.style.left=x+'px';w.style.top=y+'px';W.appendChild(w);}
 const wl=document.createElement('div');wl.className='wall';wl.style.left=G.WALLX+'px';wl.style.top=(y-30)+'px';wl.style.height='60px';W.appendChild(wl);
 const cv=document.createElement('div');cv.className='well cav';cv.style.left=G.CAVX+'px';cv.style.top=y+'px';W.appendChild(cv);});
for(let r=0;r<CY.length-1;r++){for(const x of G.JCOL){
 const l=document.createElement('div');l.className='jlink';l.style.left=x+'px';
 l.style.top=(CY[r]+30)+'px';l.style.height=(CY[r+1]-CY[r]-60)+'px';W.appendChild(l);
 const d1=document.createElement('div');d1.className='jdot';d1.style.left=x+'px';d1.style.top=(CY[r]+30)+'px';W.appendChild(d1);
 const d2=document.createElement('div');d2.className='jdot';d2.style.left=x+'px';d2.style.top=(CY[r+1]-30)+'px';W.appendChild(d2);}}
const els={};
for(const i in KIND){const e=document.createElement('div');e.className='ion '+KIND[i];
 e.textContent=KIND[i]=='yb'?'':i;els[i]=e;W.appendChild(e);}
function jkey(s){return s.length==3?null:(Math.round(s[0]*10)/10)+'|'+s[1];}
const CYs=CY;
function rowx(s,r){if(s.length==2&&typeof s[0]=='number'&&s[1]==r)return s[0];
 if(s.length==3&&s[2]==CYs[r])return s[1];return null;}
function verify(fi){return VERDICTS[fi];}  // the scheduler's per-frame verdict, precomputed
let cur=0,timer=null;
const cap=document.getElementById('cap'),ver=document.getElementById('verify'),
 sl=document.getElementById('sl'),ctr=document.getElementById('ctr');
sl.max=FR.length-1;
function show(k){cur=Math.max(0,Math.min(FR.length-1,k));const f=FR[cur];
 for(const i in f.pos){els[i].style.left=f.pos[i][0]+'px';els[i].style.top=f.pos[i][1]+'px';
  els[i].classList.toggle('hi',f.hi.includes(i));}
 cap.innerHTML=f.cap+(f.badge?'<span id="badge">'+f.badge+'</span>':'');
 const v=verify(cur);
 ver.innerHTML=v?'<span class="bad">&#9888; '+v+'</span>':'<span class="ok">&#10003; frame verified: capacities ok, junction columns clear, no free passing</span>';
 sl.value=cur;ctr.textContent=(cur+1)+' / '+FR.length;}
function play(){if(timer){clearInterval(timer);timer=null;document.getElementById('pl').textContent='Play';return;}
 document.getElementById('pl').textContent='Pause';
 timer=setInterval(()=>{if(cur>=FR.length-1){play();}else show(cur+1);},+document.getElementById('sp').value);}
document.getElementById('b0').onclick=()=>show(0);
document.getElementById('be').onclick=()=>show(FR.length-1);
document.getElementById('bp').onclick=()=>show(cur-1);
document.getElementById('bn').onclick=()=>show(cur+1);
document.getElementById('pl').onclick=play;
sl.oninput=e=>show(+e.target.value);
document.addEventListener('keydown',e=>{if(e.key=='ArrowRight')show(cur+1);if(e.key=='ArrowLeft')show(cur-1);if(e.key==' '){e.preventDefault();play();}});
show(0);
</script></body></html>"""


def verify(FR, g, kinds=None):
    """Reflects the scheduler's single legality pass. The rules live in
    qec_scheduler.frame_errors; the visualizer only supplies the pixel
    descriptor and adds no rule of its own."""
    return S.frame_errors(FR, g.rule_geom(), kinds)


def write_html(path, d, merge, rounds, distill, FR, ions, g):
    geom = {"MX": g.MX, "WELL": list(g.WELL.values()), "JCOL": list(g.JCOL.values()),
            "HOLDX": g.HOLDX if distill else [], "SWX": g.SWX,
            "PARKX": g.PARKX if merge else [], "SPX": g.SPX, "WALLX": g.WALLX,
            "CAVX": g.CAVX, "YBX": g.YBX, "CY": {str(r): g.CY[r] for r in range(d)}, "P": P, "D": d}
    capmap = {str(round(k, 1)): v for k, v in g.CAP.items()}
    verdicts = [None] * len(FR)
    for _fi, _msg in S.frame_errors(FR, g.rule_geom(), ions):
        if verdicts[_fi] is None:
            verdicts[_fi] = _msg
    if distill:
        title = f"Distance-{d} merge with double selection on the chip geometry, {rounds} rounds"
        sub = ("The distilled lane of thesis Section 4.3.4 on the full cell layout: two halves ferried to the "
               "gate-end hold wells during the round, the spent survivor recycled at the read to catch the third, "
               "double selection at the boundary. Op order from qec_distill.py.")
    elif merge:
        title = f"Distance-{d} remote merge on the chip geometry, {rounds} rounds"
        sub = ("Two full merge rounds on the cell layout of thesis Chapter 4: memory homes, gate strip, "
               "hold and swap wells, SPAM sites, wall, and cavities. Junctions are gate-zone only.")
    else:
        title = f"Distance-{d} error-correction round on the chip geometry"
        sub = ("One local round on the cell layout of thesis Chapter 4. The band files into the gate strip, "
               "gates fire in strip wells, ancillas convoy to the SPAM sites and back.")
    html = TEMPLATE.replace("__TITLE__", title).replace("__SUB__", sub)
    html = html.replace("__GEOM__", json.dumps(geom)).replace("__FRAMES__", json.dumps(FR))
    html = html.replace("__IONS__", json.dumps(ions)).replace("__CAP__", json.dumps(capmap))
    html = html.replace("__VERDICTS__", json.dumps(verdicts))
    with open(path, "w") as fh:
        fh.write(html)


def selftest():
    """Plant one violation of each verifier rule into a clean run and require
    the verifier to refuse it. A rule without a failing test is a blind spot,
    which is how the row-entry hole survived until an eyeball caught it."""
    import copy
    FR, ions, g = build(3, True, 2, False)
    assert not verify(FR, g, ions), "clean run must verify"
    k = len(FR) // 2
    anc = [i for i in FR[k]["slots"] if i.startswith("A")]
    dat = [i for i in FR[k]["slots"] if i.startswith("d")]
    def planted(mutate, name):
        BAD = copy.deepcopy(FR)
        mutate(BAD)
        n = len(verify(BAD, g, ions))
        print(f"  planted {name:28s} -> {'caught' if n else 'MISSED'} ({n} findings)")
        assert n, name
    def teleport(B):
        for j in range(k, k + 3):
            B[j]["slots"][anc[0]] = [g.WELL[0], 2]
    planted(teleport, "cross-row teleport")
    def freepass(B):
        r0 = [i for i in dat if B[k]["slots"][i][1] == 0 and isinstance(B[k]["slots"][i][0], float)]
        a2, b2 = r0[0], r0[1]
        for j in range(k, k + 3):
            B[j]["slots"][a2], B[j]["slots"][b2] = B[j]["slots"][b2], B[j]["slots"][a2]
    planted(freepass, "pass with no shared well")
    def overcap(B):
        tgt = [g.WELL[1], 0]
        for m, i in enumerate(dat[:5]):
            B[k]["slots"][i] = list(tgt)
    planted(overcap, "well over capacity")
    def jrest(B):
        B[k]["slots"][dat[0]] = [g.JCOL[1], 0]
    planted(jrest, "rest on a junction column")
    def jshare(B):
        B[k]["slots"][anc[0]] = ["J", g.JCOL[1], g.gapy(0, 1)]
        B[k]["slots"][anc[1]] = ["J", g.JCOL[1], g.gapy(0, 1)]
    planted(jshare, "two ions in one junction")
    def jswitch(B):
        B[k]["slots"][anc[0]] = ["J", g.JCOL[0], g.gapy(0, 1)]
        B[k + 1]["slots"][anc[0]] = ["J", g.JCOL[2], g.gapy(0, 1)]
    planted(jswitch, "junction switch mid-transit")
    global CHECK_ONLY
    CHECK_ONLY = True
    FR2, ions2, g2 = build(3, True, 2, True)
    inc = len(build.last_errs)
    CHECK_ONLY = False
    FR3, ions3, g3 = build(3, True, 2, True)
    bat = len(verify(FR3, g3, ions3))
    assert inc == bat == 0, f"checker paths disagree: incremental {inc}, batch {bat}"
    print(f"  incremental and batch checkers agree ({inc} findings each)")
    gf = next(j for j, f in enumerate(FR) if f.get("op") == "gate")
    def gate3(B):
        pr = B[gf]["pairs"][0]
        w = B[gf]["slots"][pr[0]]
        extra = [i for i in dat if B[gf]["slots"][i] != w][0]
        B[gf]["slots"][extra] = list(w)
    planted(gate3, "third ion in a gating well")
    def gatemem(B):
        pr = B[gf]["pairs"][0]
        for i in pr:
            B[gf]["slots"][i] = [g.MX[2], 0]
    planted(gatemem, "gate outside the gate zone")
    rf = next(j for j, f in enumerate(FR) if f.get("op") == "read")
    def readaway(B):
        pr = B[rf]["pairs"][0]
        B[rf]["slots"][pr[0]] = [g.MX[2], 0]
    planted(readaway, "read away from SPAM or swap well")
    def commmem(B):
        B[k]["slots"]["C0B"] = [g.MX[2], 0]
    planted(commmem, "comm ion in the memory zone")
    def datapastwall(B):
        B[k]["slots"][dat[0]] = [g.CAVX, 0]
    planted(datapastwall, "code ion past the optical wall")
    print("selftest: every planted violation caught")


if __name__ == "__main__":
    import qec_scheduler as S
    import qec_distill
    if len(sys.argv) > 1 and sys.argv[1] == "selftest":
        selftest()
        sys.exit(0)
    sweep = len(sys.argv) > 1 and sys.argv[1] == "all"
    ds = list(range(3, 29, 2)) if sweep else ([int(sys.argv[1])] if len(sys.argv) > 1 else [3, 5])
    if sweep:
        CHECK_ONLY = True                            # slot-level sweep: verify every distance, write no pages
    rows = []
    for d in ds:
        nchk = S.certify(d)
        qec_distill.certify_distill(d)
        cells = {}
        for merge, rounds, tag, dst in [(False, 1, "round", False),
                                        (True, MERGE_ROUNDS, "merge_full", False),
                                        (True, MERGE_ROUNDS, "merge_distill", True)]:
            FR, ions, g = build(d, merge, rounds, dst)
            errs = build.last_errs if CHECK_ONLY else verify(FR, g, ions)
            path = os.path.join(HERE, f"qec_{tag}_sim_d{d}.html")
            if not CHECK_ONLY:
                write_html(path, d, merge, rounds, dst, FR, ions, g)
            cells[tag] = (len(parallel_steps(d, merge, rounds)), len(FR), len(errs))
            print(f"d={d} {tag:13s}: {len(FR)} frames, {len(errs)} findings"
                  + ("" if CHECK_ONLY else f" -> {os.path.basename(path)}"))
            for item in errs[:6]:
                fi, e = item if isinstance(item, tuple) else ("?", item)
                print(f"   frame {fi}: {e}")
        rows.append((d, nchk, cells["round"], cells["merge_full"], cells["merge_distill"]))
    if sweep:
        print()
        print("| d | scheduler checks | local round steps | frames | findings | 2-round merge steps | frames | findings | distilled merge frames | findings |")
        print("|--:|:----------------:|------------------:|-------:|---------:|--------------------:|-------:|---------:|-----------------------:|---------:|")
        for d, nchk, (rs, rf, ro), (ms, mf, mo), (zs, zf, zo) in rows:
            print(f"| {d} | {nchk} PASS | {rs} | {rf} | {ro} | {ms} | {mf} | {mo} | {zf} | {zo} |")

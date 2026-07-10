"""Certified d=3 merge schedule replayed on the Chapter-4 geometry, v2.
Slot model: every ion occupies exactly one rest slot (a well) or a junction
transit point. Co-residency in one well = one crystal (legal up to the well's
capacity). Order along a row may change only between ions that shared a well.
The verifier enforces capacity, junction-rest, and the no-free-passing rule."""
import json, sys
import os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from qec_scheduler import build_stabilizers, place, stab_cell, num, seam_schedule, round_ops

D = 3
STABS = build_stabilizers(D)
NAMES = {frozenset({1,2,4,5}):"A1", frozenset({2,3,5,6}):"A2",
         frozenset({4,5,7,8}):"A3", frozenset({5,6,8,9}):"A4",
         frozenset({2,3}):"B1", frozenset({7,8}):"B2",
         frozenset({1,4}):"B3", frozenset({6,9}):"B4"}
LABEL = {s: NAMES[frozenset(num(rc,D) for rc in s.data)] for s in STABS}
SEAM = seam_schedule(D); CHAINS = place(D)
ID_DATA = {(r,c): "d%d"%(r*D+c+1) for r in range(D) for c in range(D)}
def A(s): return LABEL[s]
def Dt(rc): return ID_DATA[rc]

P = 46.0
MX  = [120 + k*P for k in range(7)]
GX  = [520 + m*P for m in range(6)]
JCOL = {0: GX[0], 1: GX[2], 2: GX[4]}
WELL = {0: GX[1], 1: GX[3], 2: GX[5]}
HOLDX = [796.0, 842.0, 888.0]
SWX = 934.0
PARKX = [990.0]
SPX = [1064 + k*P for k in range(3)]
WALLX = 1210.0
CAVX = 1264.0; YBX = 1296.0
CY = {r: 110 + 190*r for r in range(D)}
IW = 9.0
def gapy(a,b): return (CY[a]+CY[b])/2
SLOTS = MX + list(WELL.values()) + HOLDX + [SWX] + PARKX + SPX + [CAVX, YBX]
CAP = {**{x:2 for x in MX}, **{x:4 for x in WELL.values()}, **{x:2 for x in HOLDX}, SWX:3,
       **{x:2 for x in PARKX}, **{x:2 for x in SPX}, CAVX:2, YBX:1}

ions, home = {}, {}
for ci, chain in enumerate(CHAINS):
    for j,(kind,item) in enumerate(chain):
        i = ID_DATA[item] if kind=="data" else LABEL[item]
        home[i] = (MX[j], ci)
for r in range(D):
    for c in range(D):
        ions[ID_DATA[(r,c)]] = "data"
for s in STABS: ions[LABEL[s]] = "X" if s.kind=="X" else "Z"
for r in range(D):
    ions["C%d"%r] = "comm" if r in SEAM else "spare"
    home["C%d"%r] = (CAVX, r)
    if r in SEAM:
        ions["C%dB"%r] = "comm"
        home["C%dB"%r] = (SWX, r) if (len(sys.argv)<2 or sys.argv[1]!="round") else (CAVX, r)
        if len(sys.argv)>1 and sys.argv[1]=="distill":
            for k in range(3):
                hid = "h%d" % (3*r+k+1)
                ions[hid] = "held"; home[hid] = (HOLDX[k], r)
    ions["Y%d"%r] = "yb"; home["Y%d"%r] = (YBX, r)

slot = {i: home[i] for i in home}          # ion -> (slot_x, row) or ("J",c,gapy)
FR = []
ROT_THROUGH_SWAP = []

def groups(row):
    g = {}
    for i,s in slot.items():
        if isinstance(s, tuple) and len(s)==2 and s[1]==row and isinstance(s[0], float):
            g.setdefault(s[0], []).append(i)
    return g

def render():
    posr = {}
    for i,s in slot.items():
        if isinstance(s, tuple) and len(s)==3:
            posr[i] = [s[1], s[2]]
        else:
            x,row = s
            mem = sorted(groups(row).get(x, [i]))
            n = len(mem); k = mem.index(i)
            off = 0 if n==1 else (k - (n-1)/2) * 2*IW
            posr[i] = [round(x+off,1), CY[row]]
    return posr

def snap(cap, hi=None, junc=None, badge=""):
    FR.append({"pos": render(), "cap": cap, "hi": hi or [], "junc": junc or [],
               "badge": badge,
               "slots": {i: (list(slot[i]) if len(slot[i])==3 else [slot[i][0], slot[i][1]]) for i in slot}})

def occ(x,row): return [i for i,s in slot.items() if len(s)==2 and s[0]==x and s[1]==row]

def hop(i, tx, row, cap, hi=None, badge=""):
    """Move ion i from its slot to slot tx on the same row, rotating through
    every occupied slot on the way (each rotation = one shared-well frame)."""
    x0 = slot[i][0]; step = 1 if tx>x0 else -1
    path = sorted([x for x in SLOTS if (x-x0)*step>0.5 and (tx-x)*step>0.5], key=lambda x:(x-x0)*step)
    for bx in path:
        if occ(bx,row):
            slot[i] = (bx,row)
            snap(cap + f" {i} merges into the occupied well on its path and the crystal rotates it through.", hi=hi or [i], badge=badge or "rotate-through")
            if bx == SWX: ROT_THROUGH_SWAP.append((len(FR)-1, i))
    slot[i] = (tx,row)
    snap(cap, hi=hi or [i], badge=badge)

def multi(assign, cap, hi=None, badge=""):
    """One frame: several ions re-slot at once (used only for order-preserving
    convoys, checked by the verifier)."""
    for i,(x,row) in assign.items(): slot[i] = (x,row)
    snap(cap, hi=hi or list(assign), badge=badge)

BANDW = {}
def band_in():
    global BANDW
    BANDW = {}
    assign = {}
    for row in range(D):
        chain = sorted([i for i in slot if len(slot[i])==2 and slot[i][1]==row and slot[i][0] in MX],
                       key=lambda i: slot[i][0])
        ch = [c for c in chain if ions[c] in ("data","X","Z")]
        gs = [ch[k:k+2] for k in range(0,len(ch),2)]
        for w,g in enumerate(gs):
            for i in g:
                assign[i] = (WELL[min(w,2)], row); BANDW[i] = WELL[min(w,2)]
    multi(assign, "Gate feed: each row's band files into the gate strip in order, two ions to a well. One entry, as priced.", badge="band enters strip")

def band_out():
    assign = {i: home[i] for i in BANDW if len(slot[i])==2 and slot[i][0] in list(WELL.values())}
    multi(assign, "The band files back to its memory homes, one return, order preserved.", badge="band returns")

def to_well(i, wx, row, cap, hi=None):
    if slot[i] != (wx,row):
        hop(i, wx, row, cap, hi=hi)

MODE = sys.argv[1] if len(sys.argv)>1 else "merge"
_ops = round_ops(D, MODE in ("merge","distill"), 2)
net = {l:"C%d"%l for l in SEAM}; car = {l:"C%dB"%l for l in SEAM}
surv = {}
def C(l): return surv[l] if MODE=="distill" and l in surv else car[l]

PEND = {}

def free_hold_bodies(l):
    return sorted([i for i in ions if ions[i]=="held" and slot[i][1]==l and slot[i][0] in HOLDX],
                  key=lambda i: slot[i][0])

def ferry_half(l, k, hh, tag):
    snap(f"{tag} lane {l}: the networker heralds raw pair {k} at the cavity.", hi=[net[l]], badge="herald")
    hop(net[l], SWX, l, f"The networker carries the fresh half to the swap well.")
    snap(f"Crystal rotation at the swap well: the carrier {car[l]} takes the half.", hi=[net[l], car[l]])
    hop(net[l], CAVX, l, "The networker returns to the cavity and keeps attempting.")
    hop(car[l], slot[hh][0], l, f"The carrier parks half {k} beside {hh} in a gate-end hold well, a one-hop ferry.")
    hop(car[l], SWX, l, "The carrier returns to the swap well.")

def distill_batch(l, order, tag):
    sv, c1, c2 = order
    hop(sv, WELL[1], l, f"{tag} lane {l}: {sv}, the first catch, takes the survivor's middle gate well.")
    hop(c1, WELL[0], l, f"{c1} takes the left gate well.")
    hop(c2, WELL[2], l, f"{c2} takes the right gate well.")
    slot[c1] = (WELL[1], l)
    snap(f"Bilateral CNOT one: {c1} merges with the survivor and the gate fires.", hi=[c1, sv])
    slot[c1] = (WELL[0], l)
    snap(f"{c1} splits back out.", hi=[c1])
    hop(c1, SPX[0], l, f"{c1} hops into the SPAM zone and is read in Z, the bit-flip check.", badge="read Z")
    slot[c2] = (WELL[1], l)
    snap(f"Bilateral CNOT two: {c2} merges with the survivor while {c1} is still being read.", hi=[c2, sv])
    slot[c2] = (WELL[2], l)
    snap(f"{c2} splits back out and takes its basis rotation.", hi=[c2])
    hop(c2, SPX[1], l, f"{c2} is read in X, the phase-flip check.", badge="read X")
    snap(f"Both modules' checks agree: keep. {sv} holds the purified pair.", hi=[sv], badge="keep")
    hop(c1, HOLDX[1], l, f"{c1} resets and returns to the hold pool.")
    hop(c2, HOLDX[2], l, f"{c2} resets and returns to the hold pool.")
    hop(sv, SWX, l, f"{sv} moves to the swap well, purified pair ready for the seam.")
    surv[l] = sv
in_band = False
FILEOUT = {}
for op in _ops:
    v = op[0]
    if v == "prep":
        snap("Rest. Memory homes hold the code row" + (", the carrier waits at each active swap well, the networker and Yb sit at each cavity." if MODE!="round" else ". The comm ions rest at their cavities; a local round asks nothing of them."), badge="chapter-4 geometry")
        if MODE == "distill":
            snap("Warm-up: before round 1 each active lane distills its first batch.", badge="tier 3")
            for l in sorted(SEAM):
                pool = free_hold_bodies(l)
                for k, hh in enumerate(pool):
                    ferry_half(l, k, hh, "Warm-up")
                distill_batch(l, pool, "Warm-up")
    elif v == "round":
        snap(f"Round {op[1]+1} of 2.", badge=f"round {op[1]+1}")
    elif v == "park":
        for s, cl, well in op[1]:
            hop(A(s), MX[6], cl, f"{A(s)} idles this merge and moves to the row's end well.")
            slot[A(s)] = ("J", JCOL[D-1], gapy(min(cl,D-1), min(cl,D-1)+1) if cl < D-1 else CY[cl])
            snap(f"{A(s)} lifts into the boundary junction and descends toward the bottom cell.", hi=[A(s)], junc=[[D-1,k] for k in range(cl,D-1)], badge="in transit")
            slot[A(s)] = (PARKX[well], D-1)
            snap(f"{A(s)} drops into the bottom cell's spare park well.", hi=[A(s)])
    elif v == "unpark":
        for s, cl, well in reversed(op[1]):
            slot[A(s)] = ("J", JCOL[D-1], gapy(cl, cl+1) if cl < D-1 else CY[cl])
            snap(f"{A(s)} lifts out of its park well into the boundary junction.", hi=[A(s)], junc=[[D-1,k] for k in range(cl,D-1)], badge="in transit")
            slot[A(s)] = (MX[6], cl)
            snap(f"{A(s)} lands back on its home row's end well.", hi=[A(s)])
            hop(A(s), home[A(s)][0], cl, f"{A(s)} walks back to its memory home.")
    elif v in ("inrow","swap","xlift","xgate","xlower","xdrop","comm_out","comm_lift","comm_gate","comm_lower","comm_back","comm_arrive"):
        if not in_band:
            band_in(); in_band = True
        if v == "inrow":
            L = op[1]+1
            for s, rc in op[2]:
                if slot[Dt(rc)] != slot[A(s)]:
                    to_well(Dt(rc), slot[A(s)][0], slot[A(s)][1], f"Step {L}: {Dt(rc)} steps into {A(s)}'s gate well.")
            snap(f"Step {L}: in-row gates fire, each pair one crystal in a strip gate well.",
                 hi=[A(s) for s,_ in op[2]])
        elif v == "swap":
            L = op[1]+1
            for s, rc in op[2]:
                a, dd = A(s), Dt(rc)
                sa, sd = slot[a], slot[dd]
                slot[a] = sd
                snap(f"Step {L}: {a} merges into {dd}'s well to pass it.", hi=[a])
                slot[a], slot[dd] = sd, sa
                snap(f"Step {L}: the crystal rotates and splits; {a} and {dd} have exchanged wells.", hi=[a])
        elif v == "xlift":
            L = op[1]+1
            for s,c,ac,tr in op[2]:
                slot[A(s)] = ("J", JCOL[c], gapy(min(ac,tr),max(ac,tr)))
            snap(f"Step {L}: cross-row ancillas lift into their junction mouths at c + 1/4.",
                 hi=[A(s) for s,*_ in op[2]], junc=[[c,min(ac,tr)] for s,c,ac,tr in op[2]], badge="in transit")
        elif v == "xgate":
            L = op[1]+1
            for s, rc in op[2]:
                tw = slot[Dt(rc)]
                slot[A(s)] = tw
            snap(f"Step {L}: each drops into the neighbor row's strip well and gates its data there.",
                 hi=[A(s) for s,_ in op[2]])
        elif v == "xlower":
            L = op[1]+1
            for s,c,ac,tr in op[2]:
                slot[A(s)] = ("J", JCOL[c], gapy(min(ac,tr),max(ac,tr)))
            snap(f"Step {L}: they lift back into their junction mouths.",
                 hi=[A(s) for s,*_ in op[2]], junc=[[c,min(ac,tr)] for s,c,ac,tr in op[2]], badge="in transit")
        elif v == "xdrop":
            L = op[1]+1
            for s,c,ac,swapped in op[2]:
                slot[A(s)] = (BANDW.get(A(s), WELL[2]), ac)
            snap(f"Step {L}: they drop back onto their own strip wells.", hi=[A(s) for s,*_ in op[2]])
        elif v == "comm_out":
            L = op[1]+1
            for l in op[2]:
                hop(C(l), WELL[2], l, f"Step {L}: the seam ion steps from the swap well into the strip, toward its boundary data.", badge="seam excursion")
        elif v == "comm_lift":
            L = op[1]+1
            for l in op[2]:
                slot[C(l)] = ("J", JCOL[D-1], gapy(l,l+1))
            snap(f"Step {L}: the seam ions lift into the boundary junction mouths.", hi=[C(l) for l in op[2]],
                 junc=[[D-1,l] for l in op[2]], badge="in transit")
        elif v == "comm_gate":
            L = op[1]+1
            for l, rc in op[2]:
                slot[C(l)] = slot[Dt(rc)]
            snap(f"Step {L}: the seam gate fires in the strip well holding the boundary data.",
                 hi=[C(l) for l,_ in op[2]])
        elif v == "comm_lower":
            L = op[1]+1
            for l in op[2]:
                slot[C(l)] = ("J", JCOL[D-1], gapy(l,l+1))
            snap(f"Step {L}: the seam ions lift back into the junction mouths.", hi=[C(l) for l in op[2]],
                 junc=[[D-1,l] for l in op[2]], badge="in transit")
        elif v == "comm_back":
            L = op[1]+1
            for l in op[2]:
                slot[C(l)] = (WELL[2], l)
                snap(f"Step {L}: the carrier drops back onto its own strip.", hi=[C(l)])
                hop(C(l), SWX, l, f"Step {L}: it returns to the gate-end swap well to wait for the round's read.", badge="seam excursion")
        elif v == "comm_arrive":
            for l in op[2]:
                if slot[C(l)] != (SWX,l):
                    hop(C(l), SWX, l, "The carrier settles at its swap well.")
    elif v == "readout":
        if in_band: band_out(); in_band = False
        pairs = [(A(s), Dt(rc)) for s,rc in op[1]]
        for a,dd in pairs:
            sa, sd = slot[a], slot[dd]
            slot[a] = sd
            snap(f"File-out: {a} merges into {dd}'s well.", hi=[a])
            slot[a], slot[dd] = sd, sa
            snap(f"The crystal rotates and splits: {a} has passed {dd}.", hi=[a])
    elif v == "to_spam":
        order = sorted([A(s) for s in op[1]], key=lambda i: -slot[i][0])
        FILEOUT.update({a: slot[a] for a in order})
        for a in order:
            row = slot[a][1]
            free = [x for x in SPX if not occ(x,row)]
            hop(a, free[0], row, f"{a} shuttles the corridor, across the strip, to a SPAM site.", badge="convoy")
        snap("All ancillas sit at SPAM detection sites, behind the wall from the memory data.", hi=order, badge="in SPAM")
    elif v == "syndromes":
        snap("493 nm readout fires at the SPAM sites" + (", the seam bits read at the swap wells in the same beat." if op[2] else "."),
             hi=[A(s) for s in op[1]], badge="syndromes")
    elif v == "measure":
        if MODE == "distill":
            snap("The survivors are read at the swap wells. This round's purified pairs reach module B.",
                 hi=[C(l) for l in op[1]], badge="purified pairs to B")
            for l in op[1]:
                old = surv.pop(l, None)
                if old is None:
                    continue
                free = [x for x in HOLDX if not occ(x,l)]
                hop(old, free[0], l, f"{old} resets and rejoins the hold pool, an empty body again.")
                ferry_half(l, 2, old, "Third catch,")
                distill_batch(l, PEND.get(l, [])[:2] + [old], "Boundary:")
        else:
            snap("The carriers are read at their swap wells. This round's Bell pairs reach module B.",
                 hi=[C(l) for l in op[1]], badge="pairs to B")
    elif v == "from_spam":
        order = sorted([A(s) for s in op[1]], key=lambda i: FILEOUT[A(s)][0] if A(s) in FILEOUT else slot[i][0])
        order = sorted([A(s) for s in op[1]], key=lambda a: -FILEOUT[a][0])
        for a in order:
            tx, row = FILEOUT[a]
            hop(a, tx, row, f"{a} shuttles back in from the SPAM zone to where it left the row.", badge="convoy")
    elif v == "reset":
        pairs = [(A(s), Dt(rc)) for s,rc in op[1]]
        for a,dd in pairs:
            sa, sd = slot[a], slot[dd]
            slot[a] = sd
            snap(f"Reset file-in: {a} merges back into {dd}'s well.", hi=[a])
            slot[a], slot[dd] = sd, sa
            snap(f"Rotation undone: {a} is back on the memory side of {dd}.", hi=[a])
    elif v == "reset_done":
        assign = {}
        for s in STABS:
            i = A(s)
            if slot[i][0] in MX or slot[i][0] in SPX:
                assign[i] = home[i]
        for n in range(1, D*D+1):
            i = "d%d"%n
            assign[i] = home[i]
        multi(assign, "Everyone is back in the memory home it started from. The round ends at rest.", badge="at rest")
    elif v == "herald":
        if MODE == "distill":
            for l in sorted(op[1]):
                pool = free_hold_bodies(l)[:2]
                for k, hh in enumerate(pool):
                    ferry_half(l, k, hh, "Next batch")
                PEND[l] = pool
            snap("Two of the next batch's halves are parked; the third body is still out at the seam. "
                 "Four live states and one attempter: the five-ion lane at full occupancy.",
                 badge="occupancy 5")
        else:
            snap("The networkers herald fresh Bell pairs at their cavities. Yb holds them cold.",
                 hi=[net[l] for l in op[1]], badge="herald")
    elif v == "comm_swap":
        if MODE == "distill":
            snap("Ping-pong roles hold: the networker keeps the cavity, the carrier keeps the ferry. The survivor pipeline replaces the raw handoff.", badge="distilled lane")
            continue
        for l in op[1]:
            hop(net[l], SWX, l, "Ping-pong: the heralded networker crosses the wall and the SPAM row to the swap well.", badge="handoff")
            snap("The two comm ions share the swap well and the crystal rotates: roles exchange.", hi=[net[l], car[l]])
            hop(car[l], CAVX, l, "The emptied ion returns to the cavity as the next networker.")
            net[l], car[l] = car[l], net[l]
    else:
        raise ValueError(v)

# ------------------------------ verifier ------------------------------------
def verify():
    errs = []
    for fi,f in enumerate(FR):
        g = {}
        for i,s in f["slots"].items():
            if len(s)==2 and not isinstance(s[0],str):
                g.setdefault((s[0],s[1]),[]).append(i)
            elif len(s)==3:
                pass
        for (x,row),mem in g.items():
            if x in CAP and len(mem) > CAP[x]:
                errs.append((fi, f"well x={round(x)} row {row} holds {len(mem)} > cap {CAP[x]}: {mem}"))
            if x in JCOL.values():
                errs.append((fi, f"{mem} at REST on junction column x={round(x)}"))
    for fi in range(1,len(FR)):
        p, f = FR[fi-1], FR[fi]
        share = set()
        for fr in (p,f):
            g = {}
            for i,s in fr["slots"].items():
                if len(s)==2 and not isinstance(s[0],str): g.setdefault((s[0],s[1]),[]).append(i)
            for mem in g.values():
                for a in mem:
                    for b in mem:
                        if a<b: share.add((a,b))
        for row in range(D):
            def seq(fr):
                lst = [(s[0],i) for i,s in fr["slots"].items() if len(s)==2 and not isinstance(s[0],str) and s[1]==row]
                return [i for _,i in sorted(lst)]
            sa, sb = seq(p), seq(f)
            rank = {i:k for k,i in enumerate(sb)}
            common = [i for i in sa if i in rank]
            for ii in range(len(common)):
                for jj in range(ii+1,len(common)):
                    a,b = common[ii], common[jj]
                    if rank[a] > rank[b] and (min(a,b),max(a,b)) not in share:
                        errs.append((fi, f"row {row}: {a} passed {b} with no shared well"))
    return errs

errs = verify()
print(f"frames: {len(FR)}  findings: {len(errs)}  rotations-through-swap-well: {len(ROT_THROUGH_SWAP)}")
seen = set()
for fi,e in errs:
    key = e.split(":")[0]+e[:40]
    if key in seen: continue
    seen.add(key)
    print(f"  f{fi:3d}: {e}   | {FR[fi]['cap'][:60]}")
    if len(seen)>=15: break
json.dump({"frames":FR,"ions":ions,"mode":MODE,
           "geom":{"MX":MX,"WELL":list(WELL.values()),"JCOL":list(JCOL.values()),"HOLDX":HOLDX,"SWX":SWX,
                   "PARKX":PARKX,"SPX":SPX,"WALLX":WALLX,"CAVX":CAVX,"YBX":YBX,
                   "CY":CY,"P":P,"D":D}}, open(os.path.join(HERE, f"chip_frames_{MODE}.json"),"w"))
print("json written")

import subprocess
subprocess.run([sys.executable, os.path.join(HERE, 'qec_chip_emit.py'), MODE], check=True)
print(f'qec_chip_sim_d3_{MODE}.html written')

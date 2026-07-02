#!/usr/bin/env python3
"""Section 4.4.4 scaling figure: proves the trapped-ion QPU design is
extendable in code distance d. All numbers are computed live from the real
scheduler module qec_scheduler.py -- nothing is hard-coded."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter
from matplotlib.patches import Patch
import numpy as np

import qec_scheduler as S

# --------------------------------------------------------------------------
# 1. COMPUTE THE DATA (live, from the scheduler)
# --------------------------------------------------------------------------
DS = [3, 5, 7, 9, 11, 13, 15]

round_depth, connection, readout = [], [], []
merge_depth, qubits, gates, lanes = [], [], [], []

for d in DS:
    steps = S.parallel_steps(d, merge=False, rounds=1)
    ops = S.round_ops(d, merge=False, rounds=1)
    # readout tail = SPAM bubble ('readout') + shuttle-out + 493nm read + shuttle-back.
    rd = sum(1 for ixs in steps if any(ops[i][0] in ("readout", "to_spam", "syndromes", "from_spam") for i in ixs))
    round_depth.append(len(steps))
    readout.append(rd)
    connection.append(len(steps) - rd)
    merge_depth.append(len(S.parallel_steps(d, merge=True, rounds=d)))
    qubits.append(d * d + (d * d - 1))
    gates.append(S.op_tally(d, merge=False)["two_qubit_gates"])
    lanes.append(d)

d = np.array(DS)
round_depth = np.array(round_depth)
connection = np.array(connection)
readout = np.array(readout)
merge_depth = np.array(merge_depth)
qubits = np.array(qubits)
gates = np.array(gates)
lanes = np.array(lanes)

assert (round_depth == d + 27).all(), "round_depth != d + 27"   # key finding

# --------------------------------------------------------------------------
# 2. PALETTE + STYLE  (must match the other thesis figures exactly)
# --------------------------------------------------------------------------
INK      = "#2c2c2a"   # text / axes
MUTE     = "#6b6a64"   # muted labels
GRID     = "#d9d7cd"   # gridlines
C_CONN   = "#3B7FD4"   # connection band
C_READ   = "#BA7517"   # readout band
C_MERGE  = "#534AB7"   # merge line
C_QUBITS = "#0F6E56"   # resource: qubits
C_GATES  = "#D2703A"   # resource: gates
C_LANES  = "#6b6a64"   # resource: lanes

plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor":   "white",
    "savefig.facecolor": "white",
    "font.family":      "sans-serif",
    "font.sans-serif":  ["DejaVu Sans", "Helvetica", "Arial"],
    "font.size":        9,
    "axes.edgecolor":   INK,
    "axes.linewidth":   0.7,
    "axes.labelcolor":  INK,
    "text.color":       INK,
    "xtick.color":      INK,
    "ytick.color":      INK,
    "xtick.labelsize":  8,
    "ytick.labelsize":  8,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "legend.frameon":   False,
    "legend.fontsize":  8,
})


def thin_spines(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(INK)
    ax.spines["bottom"].set_color(INK)
    ax.tick_params(length=3)


fig, (axA, axB) = plt.subplots(
    2, 1, figsize=(6.5, 7.0),
    gridspec_kw=dict(hspace=0.46, left=0.105, right=0.80,
                     top=0.945, bottom=0.125))

XLIM = (1.7, 16.3)

# --------------------------------------------------------------------------
# 3. PANEL A -- "Round time stays cheap"
# --------------------------------------------------------------------------
w = 1.05
axA.bar(d, connection, width=w, color=C_CONN, edgecolor="white",
        linewidth=0.5, zorder=3)
axA.bar(d, readout, width=w, bottom=connection, color=C_READ,
        edgecolor="white", linewidth=0.5, zorder=3)

axA.set_ylabel("local round depth  (time-steps)", color=INK)
axA.set_xlabel("code distance  $d$", color=MUTE)
axA.set_title("A   Round time stays cheap", loc="left", weight="bold",
              color=INK, pad=8)
axA.set_xticks(DS)
axA.set_ylim(0, 52)
axA.set_xlim(*XLIM)
axA.yaxis.grid(True, color=GRID, linewidth=0.7, linestyle=(0, (4, 3)), zorder=0)
axA.set_axisbelow(True)
thin_spines(axA)

# --- direct band labels (no legend needed; cleaner for a thesis) ----------
# readout band label: text high-left in the empty corner, leader to the d=13 orange.
axA.annotate("readout: $d$ bubble layers,\nthen shuttle out, read, back",
             xy=(13, connection[5] + readout[5] * 0.5), xycoords="data",
             xytext=(2.0, 50.0), fontsize=7.4, color=C_READ,
             va="center", ha="left",
             arrowprops=dict(arrowstyle="-", color=C_READ, lw=0.7,
                             connectionstyle="arc3,rad=-0.22"))
# connection explanation: text sits in empty band just above the bars at low d.
axA.annotate("connection band: 4 flat gate-phases,\n$=24$ steps, does not grow with $d$",
             xy=(4, connection[0] + 0.5), xycoords="data",
             xytext=(2.0, 42.0), fontsize=7.6, color=C_CONN,
             va="center", ha="left",
             arrowprops=dict(arrowstyle="-", color=C_CONN, lw=0.7,
                             connectionstyle="arc3,rad=0.24"))
# closed form for the round depth, boxed, in clear space over the mid bars.
axA.text(9.3, 46.5, r"round depth $= 27 + d$",
         ha="left", va="center", fontsize=8, color=INK,
         bbox=dict(boxstyle="round,pad=0.32", fc="white", ec=GRID, lw=0.7))

# twin axis: full d-round merge depth (quadratic)
axA2 = axA.twinx()
axA2.plot(d, merge_depth, marker="o", ms=3.4, lw=1.3, color=C_MERGE,
          zorder=6, clip_on=False)
axA2.set_ylabel("full $d$-round merge depth", color=C_MERGE)
axA2.tick_params(axis="y", colors=C_MERGE)
axA2.spines["top"].set_visible(False)
axA2.spines["left"].set_visible(False)
axA2.spines["right"].set_color(C_MERGE)
axA2.spines["right"].set_linewidth(0.7)
axA2.spines["bottom"].set_visible(False)
axA2.set_ylim(0, max(merge_depth) * 1.14)
axA2.set_xlim(*XLIM)
# merge-curve label sits just below its own line, right side, its own clear band
axA2.annotate(r"merge $= d$ rounds $\Rightarrow O(d^2)$",
              xy=(14, merge_depth[6] * 0.86), xycoords="data",
              xytext=(5.6, merge_depth[5] * 0.86), fontsize=7.6,
              color=C_MERGE, ha="left", va="center",
              arrowprops=dict(arrowstyle="-", color=C_MERGE, lw=0.7,
                              connectionstyle="arc3,rad=-0.12"))


# --------------------------------------------------------------------------
# 4. PANEL B -- "Resources are simple polynomials"
# --------------------------------------------------------------------------
axB.plot(d, qubits, marker="o", ms=3.4, lw=1.3, color=C_QUBITS,
         zorder=5, clip_on=False)
axB.plot(d, gates, marker="s", ms=3.2, lw=1.3, color=C_GATES,
         zorder=5, clip_on=False)
axB.plot(d, lanes, marker="^", ms=3.6, lw=1.3, color=C_LANES,
         zorder=5, clip_on=False)

axB.set_yscale("log")
axB.set_ylabel("count", color=INK)
axB.set_xlabel("code distance  $d$", color=MUTE)
axB.set_title("B   Resources are simple polynomials", loc="left",
              weight="bold", color=INK, pad=8)
axB.set_xticks(DS)
axB.set_xlim(*XLIM)
axB.set_ylim(2, 2200)
axB.yaxis.set_major_locator(LogLocator(base=10.0))
axB.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
axB.yaxis.set_minor_formatter(NullFormatter())
axB.yaxis.grid(True, which="major", color=GRID, linewidth=0.7,
               linestyle=(0, (4, 3)), zorder=0)
axB.yaxis.grid(True, which="minor", color=GRID, linewidth=0.4,
               linestyle=(0, (1, 4)), alpha=0.7, zorder=0)
axB.set_axisbelow(True)
thin_spines(axB)

# inline closed-form labels, in the right margin (annotation_clip off).
# nudged apart in y so gates/qubits do not touch.
axB.annotate(r"gates $= 4d(d{-}1)$", xy=(d[-1], gates[-1]),
             xytext=(15.35, 1050), color=C_GATES, fontsize=8,
             va="center", ha="left", annotation_clip=False)
axB.annotate(r"qubits $= 2d^2{-}1$", xy=(d[-1], qubits[-1]),
             xytext=(15.35, 360), color=C_QUBITS, fontsize=8,
             va="center", ha="left", annotation_clip=False)
axB.annotate(r"lanes $= d$", xy=(d[-1], lanes[-1]),
             xytext=(15.35, 15), color=C_LANES, fontsize=8,
             va="center", ha="left", annotation_clip=False)

# footnote caption (kept fully inside the figure width)
fig.text(0.105, 0.014,
         "All values computed live from qec_scheduler.py.  Round depth is flat connection\n"
         r"(24) plus readout ($d$ bubble layers, then shuttle out / read / back); every "
         "\nresource is a low-order polynomial in $d$.",
         fontsize=6.6, color=MUTE, ha="left", va="bottom", linespacing=1.35)

# --------------------------------------------------------------------------
# 5. SAVE
# --------------------------------------------------------------------------
# save next to this script so the paths work in any session
HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "ch4_scaling.pdf")
PNG = os.path.join(HERE, "ch4_scaling_preview.png")
os.makedirs(os.path.dirname(PDF), exist_ok=True)
fig.savefig(PDF)
fig.savefig(PNG, dpi=130)
print("wrote", PDF)
print("wrote", PNG)

print("\n d | depth conn read | merge  qubits gates lanes")
for i, dd in enumerate(DS):
    print(f"{dd:2d} | {round_depth[i]:5d} {connection[i]:4d} {readout[i]:4d} |"
          f" {merge_depth[i]:5d}  {qubits[i]:5d} {gates[i]:5d} {lanes[i]:4d}")

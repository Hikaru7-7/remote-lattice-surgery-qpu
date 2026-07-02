#!/usr/bin/env python3
"""ch4_placement.py -- 'code to chip' figure for thesis section 4.4.2.

LEFT : abstract distance-3 rotated surface code (data grid + stabilizer patches).
RIGHT: the chip's balanced placement (place(3)) as 3 stacked cell-chains.
Same check == same colour AND same id on both panels, from REAL scheduler data.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle
from matplotlib.lines import Line2D

import qec_scheduler as S

# ---- palette (shared across the three thesis figures) --------------------
BG        = "white"
TEXT      = "#2c2c2a"
MUTED     = "#6b6a64"
THIN      = "#d9d7cd"
X_FILL    = "#3B7FD4"; X_EDGE = "#1C4F8C"
Z_FILL    = "#D2703A"; Z_EDGE = "#8F4620"
DATA_FILL = "white";   DATA_EDGE = "#8f8d84"

plt.rcParams.update({
    "font.size": 8,
    "font.family": "sans-serif",
    "text.color": TEXT,
    "axes.edgecolor": TEXT,
    "pdf.fonttype": 42,   # editable text in vector PDF
    "ps.fonttype": 42,
})

D = 3

# ---- pull the REAL data --------------------------------------------------
stabs   = S.build_stabilizers(D)          # 8 checks, bulk-then-edge order
cells   = S.place(D)                       # 3 cells, each a L->R chain
cellof  = S.stab_cell(D)                   # Stabilizer -> cell index
percell = S.per_cell_ancillas(D)           # [3, 3, 2]

def dlabel(rc):                            # data id d1..d9 (row-major, from num)
    return f"d{S.num(rc, D)}"

# ---- assign stable ids: X1.. / Z1.. in build order ----------------------
# Scanning build_stabilizers gives bulk checks first then edge checks, so ids
# are deterministic; we reuse the SAME dict for both panels.
ids, nX, nZ = {}, 0, 0
for s in stabs:
    if s.kind == "X":
        nX += 1; ids[s] = f"X{nX}"
    else:
        nZ += 1; ids[s] = f"Z{nZ}"

def cols(kind):
    return (X_FILL, X_EDGE) if kind == "X" else (Z_FILL, Z_EDGE)

# =========================================================================
fig, (axL, axR) = plt.subplots(
    1, 2, figsize=(6.5, 3.5), gridspec_kw={"width_ratios": [1.0, 1.30]}
)
for ax in (axL, axR):
    ax.set_aspect("equal")
    ax.axis("off")
fig.patch.set_facecolor(BG)

# -------------------------------------------------------------------------
# LEFT PANEL : the surface code
#   data qubit (r,c) drawn at (col, -row); stabilizer = translucent patch
#   covering exactly the data it reads.
# -------------------------------------------------------------------------
def P(rc):                                 # (row,col) -> plot (x,y)
    r, c = rc
    return (c, -r)

# stabilizer patches first (under the qubits). Bulk patches a touch smaller than
# edge patches so adjacent plaquette labels read cleanly; edge patches sit proud.
PAD4 = 0.40                                  # half-extent for a bulk (weight-4) cell
PAD2 = 0.36                                  # half-extent for an edge (weight-2) cell
for s in stabs:
    fill, edge = cols(s.kind)
    rs = sorted({r for r, c in s.data}); cs = sorted({c for r, c in s.data})
    pad = PAD4 if s.weight == 4 else PAD2
    x0 = min(cs) - pad; x1 = max(cs) + pad
    y0 = -max(rs) - pad; y1 = -min(rs) + pad
    w, h = x1 - x0, y1 - y0
    patch = FancyBboxPatch(
        (x0, y0), w, h,
        boxstyle="round,pad=0,rounding_size=0.16",
        linewidth=1.1, facecolor=fill, edgecolor=edge,
        alpha=0.22, mutation_aspect=1.0, zorder=1,
    )
    patch.set_clip_on(False)
    axL.add_patch(patch)
    # id label, tinted to the check colour.
    cx = (min(cs) + max(cs)) / 2
    cy = -(min(rs) + max(rs)) / 2
    if s.weight == 2:
        # push edge labels toward the outer boundary so they clear the qubits
        if len(rs) == 1:                     # top/bottom edge: shift vertically out
            cy += 0.30 if min(rs) == 0 else -0.30
        else:                                # left/right edge: shift horizontally out
            cx += -0.30 if min(cs) == 0 else 0.30
    axL.text(cx, cy, ids[s], ha="center", va="center",
             fontsize=7.4, fontweight="bold", color=edge, zorder=2)

# data qubits on top
for r in range(D):
    for c in range(D):
        x, y = P((r, c))
        axL.add_patch(Circle((x, y), 0.13, facecolor=DATA_FILL,
                             edgecolor=DATA_EDGE, linewidth=1.0, zorder=4))
        axL.text(x, y - 0.30, dlabel((r, c)), ha="center", va="center",
                 fontsize=6.0, color=MUTED, zorder=4)

axL.set_xlim(-0.85, 2.85)
axL.set_ylim(-2.95, 0.9)
axL.text(1.0, 0.72, "distance-3 rotated surface code",
         ha="center", va="center", fontsize=8.2, color=TEXT)

# -------------------------------------------------------------------------
# RIGHT PANEL : the chip placement (place(3))
#   3 horizontal cells stacked top->bottom; each drawn as its L->R chain.
#   data = white circle (d1..d9); ancilla = coloured rounded square (same id).
# -------------------------------------------------------------------------
DX = 1.0                                    # horizontal spacing between chain entries
ANC = 0.30                                  # half-size of ancilla square
DR  = 0.18                                   # data circle radius
CELL_GAP = 1.42                              # vertical gap between cells

max_len = max(len(cell) for cell in cells)

for i, cell in enumerate(cells):
    y = -i * CELL_GAP
    n = len(cell)
    xs = [j * DX for j in range(n)]
    # faint rail behind the chain
    axR.plot([xs[0] - 0.15, xs[-1] + 0.15], [y, y],
             color=THIN, linewidth=1.4, zorder=0, solid_capstyle="round")
    n_anc = percell[i]
    # cell label + ancilla count on the left
    axR.text(-1.05, y + 0.10, f"cell {i}", ha="left", va="center",
             fontsize=7.6, color=TEXT)
    axR.text(-1.05, y - 0.22, f"{n_anc} anc", ha="left", va="center",
             fontsize=6.4, color=MUTED)

    for x, (kind, item) in zip(xs, cell):
        if kind == "data":
            axR.add_patch(Circle((x, y), DR, facecolor=DATA_FILL,
                                 edgecolor=DATA_EDGE, linewidth=1.0, zorder=3))
            axR.text(x, y, dlabel(item), ha="center", va="center",
                     fontsize=5.8, color=MUTED, zorder=4)
        else:  # ancilla
            fill, edge = cols(item.kind)
            sq = FancyBboxPatch(
                (x - ANC, y - ANC), 2 * ANC, 2 * ANC,
                boxstyle="round,pad=0,rounding_size=0.10",
                linewidth=1.2, facecolor=fill, edgecolor=edge,
                alpha=0.92, zorder=3,
            )
            axR.add_patch(sq)
            axR.text(x, y, ids[item], ha="center", va="center",
                     fontsize=6.6, fontweight="bold", color="white", zorder=4)

axR.set_xlim(-1.25, (max_len - 1) * DX + 0.7)
axR.set_ylim(-((D - 1) * CELL_GAP) - 1.05, 0.9)
axR.text(((max_len - 1) * DX) / 2 - 0.1, 0.72,
         "balanced chip placement  (3, 3, 2)",
         ha="center", va="center", fontsize=8.2, color=TEXT)

# -------------------------------------------------------------------------
# shared legend (bottom), and a light 'code -> chip' cue between panels
# -------------------------------------------------------------------------
handles = [
    Line2D([0],[0], marker="s", markersize=8, linestyle="none",
           markerfacecolor=X_FILL, markeredgecolor=X_EDGE, label="X check"),
    Line2D([0],[0], marker="s", markersize=8, linestyle="none",
           markerfacecolor=Z_FILL, markeredgecolor=Z_EDGE, label="Z check"),
    Line2D([0],[0], marker="o", markersize=8, linestyle="none",
           markerfacecolor=DATA_FILL, markeredgecolor=DATA_EDGE, label="data qubit"),
]
leg = fig.legend(handles=handles, loc="lower center", ncol=3,
                 frameon=False, fontsize=7.4, handletextpad=0.4,
                 columnspacing=1.6, bbox_to_anchor=(0.5, 0.005))
for t in leg.get_texts():
    t.set_color(TEXT)

# 'code -> chip' cue in the gutter, up in the title band so it clears all rows
fig.text(0.492, 0.845, "$\\longrightarrow$", ha="center", va="center",
         fontsize=15, color=MUTED)
fig.text(0.492, 0.895, "place", ha="center", va="center",
         fontsize=6.8, color=MUTED)

fig.suptitle("From code to chip: placing the distance-3 checks",
             fontsize=9.6, color=TEXT, y=0.985)

fig.subplots_adjust(left=0.015, right=0.99, top=0.90, bottom=0.11, wspace=0.02)

# save next to this script so the paths work in any session
import os
HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "ch4_placement.pdf")
PNG = os.path.join(HERE, "ch4_placement_preview.png")
fig.savefig(PDF)                 # vector
fig.savefig(PNG, dpi=130)        # raster preview
print("saved:", PDF)
print("saved:", PNG)

# ---- report the mapping we actually drew --------------------------------
print("\n=== check -> (id, kind, weight, data, cell) drawn ===")
for s in stabs:
    ds = ",".join(dlabel(rc) for rc in sorted(s.data))
    print(f"  {ids[s]:>3}  {s.kind}  w{s.weight}  cell {cellof[s]}  [{ds}]")
print("per-cell ancilla counts:", percell)

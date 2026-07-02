#!/usr/bin/env python3
"""Section 4.4.3 "one merge round" figure.

Rebuilt as VECTOR output that matches the project's HTML visualizer EXACTLY.
The visualizer (qec_visualizer.py) draws each frame as SVG in its JS render();
this script ports that SAME SVG drawing to Python, composes 4 frames into a
filmstrip, and exports a vector PDF (via pycairo) plus a PNG preview.

Frame geometry is NOT remapped -- we use the module's own coordinates so the
look is pixel-faithful to the interactive visualizer.
"""
import os
import math
import cairo
import qec_visualizer as V

# ---- EXACT palette (light theme, as the HTML :root defaults) ---------------
BG      = "#faf9f5"
PANEL   = "#ffffff"
INK     = "#2c2c2a"
MUT     = "#6b6a64"
LINE    = "#d9d7cd"
C_X     = "#3B7FD4"
C_XD    = "#1C4F8C"
C_Z     = "#D2703A"
C_ZD    = "#8F4620"
PURPLE  = "#534AB7"
TEAL    = "#0F6E56"
AMBER   = "#BA7517"
COMM_DK = "#0a3a2e"

VARS = {
    "var(--bg)": BG, "var(--panel)": PANEL, "var(--ink)": INK, "var(--mut)": MUT,
    "var(--line)": LINE, "var(--x)": C_X, "var(--xd)": C_XD, "var(--z)": C_Z,
    "var(--zd)": C_ZD, "var(--purple)": PURPLE, "var(--teal)": TEAL,
    "var(--amber)": AMBER,
}


def resolve(c):
    return VARS.get(c, c)


def hex_rgb(c):
    c = resolve(c).lstrip("#")
    if len(c) == 3:
        c = "".join(ch * 2 for ch in c)
    return tuple(int(c[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


class Canvas:
    def __init__(self, ctx):
        self.ctx = ctx

    def push(self):
        self.ctx.save()

    def pop(self):
        self.ctx.restore()

    def translate(self, x, y):
        self.ctx.translate(x, y)

    def _set_source(self, color, opacity=1.0):
        r, g, b = hex_rgb(color)
        self.ctx.set_source_rgba(r, g, b, opacity)

    def rounded_rect(self, x, y, w, h, rx):
        c = self.ctx
        rx = min(rx, w / 2, h / 2)
        c.new_sub_path()
        c.arc(x + w - rx, y + rx, rx, -math.pi / 2, 0)
        c.arc(x + w - rx, y + h - rx, rx, 0, math.pi / 2)
        c.arc(x + rx, y + h - rx, rx, math.pi / 2, math.pi)
        c.arc(x + rx, y + rx, rx, math.pi, 3 * math.pi / 2)
        c.close_path()

    def rect(self, x, y, w, h, rx=0, fill=None, fill_opacity=1.0,
             stroke=None, stroke_width=1.0, opacity=1.0, dash=None):
        c = self.ctx
        if rx > 0:
            self.rounded_rect(x, y, w, h, rx)
        else:
            c.rectangle(x, y, w, h)
        if fill is not None:
            self._set_source(fill, fill_opacity * opacity)
            c.fill_preserve()
        if stroke is not None:
            self._set_source(stroke, opacity)
            c.set_line_width(stroke_width)
            c.set_dash(dash or [], 0)
            c.stroke()
        else:
            c.new_path()
        c.set_dash([], 0)

    def line(self, x1, y1, x2, y2, stroke, stroke_width=1.0, opacity=1.0,
             cap="butt", dash=None):
        c = self.ctx
        c.move_to(x1, y1)
        c.line_to(x2, y2)
        self._set_source(stroke, opacity)
        c.set_line_width(stroke_width)
        c.set_line_cap({"butt": cairo.LINE_CAP_BUTT,
                        "round": cairo.LINE_CAP_ROUND}[cap])
        c.set_dash(dash or [], 0)
        c.stroke()
        c.set_dash([], 0)
        c.set_line_cap(cairo.LINE_CAP_BUTT)

    def circle(self, cx, cy, r, fill=None, stroke=None, stroke_width=1.0,
               opacity=1.0):
        c = self.ctx
        c.arc(cx, cy, r, 0, 2 * math.pi)
        if fill is not None:
            self._set_source(fill, opacity)
            c.fill_preserve()
        if stroke is not None:
            self._set_source(stroke, opacity)
            c.set_line_width(stroke_width)
            c.stroke()
        else:
            c.new_path()

    def arc_path(self, cx, cy, r, a1, a2, stroke, stroke_width=1.0,
                 opacity=1.0):
        c = self.ctx
        c.new_sub_path()
        c.arc(cx, cy, r, a1, a2)
        self._set_source(stroke, opacity)
        c.set_line_width(stroke_width)
        c.stroke()

    def text(self, x, y, s, size, fill, anchor="start", weight="normal"):
        c = self.ctx
        c.select_font_face(
            "DejaVu Sans", cairo.FONT_SLANT_NORMAL,
            cairo.FONT_WEIGHT_BOLD if weight in ("bold", "600", 600)
            else cairo.FONT_WEIGHT_NORMAL)
        c.set_font_size(size)
        ext = c.text_extents(s)
        if anchor == "middle":
            xx = x - ext.width / 2 - ext.x_bearing
        elif anchor == "end":
            xx = x - ext.width - ext.x_bearing
        else:
            xx = x - ext.x_bearing
        self._set_source(fill, 1.0)
        c.move_to(xx, y)
        c.show_text(s)
        c.new_path()


def draw_frame(cv, DATA, frame):
    W = DATA["xhi"]
    celly = DATA["celly"]
    xlo = DATA["xlo"]
    xif = DATA["xif"]
    ions = DATA["ions"]
    pos = frame["pos"]
    hi = set(frame["hi"])
    junc_active = set(tuple(j) for j in frame["junc"])
    merged = frame["merged"]

    # 1. cell zones
    for ci, y in enumerate(celly):
        cv.rect(xlo, y - 24, W - xlo - 20, 48, rx=10, fill=None,
                stroke=LINE, stroke_width=1.0, opacity=.55, dash=[5, 4])
        cv.text(xlo + 5, y - 31, "cell %d" % ci, 11, MUT)

    # 2. interface
    if xif:
        cv.text(xif - 6, celly[0] - 31, "interface (cavities)", 11, MUT)
        for y in celly:
            cv.arc_path(xif + 30, y, 17, -math.pi / 2, math.pi / 2,
                        stroke=TEAL, stroke_width=2, opacity=.55)

    # 3. wells
    for wx, wy in DATA["wells"]:
        cv.rect(wx - 22, wy - 19, 44, 38, rx=10, fill=PANEL, stroke=LINE,
                stroke_width=1.2)

    # 4. junctions
    for j in DATA["junctions"]:
        c, b, jx, y1, y2 = j["c"], j["b"], j["x"], j["y1"], j["y2"]
        active = (c, b) in junc_active
        if active:
            cv.line(jx, y1 + 24, jx, y2 - 24, stroke=AMBER, stroke_width=4,
                    opacity=1.0, cap="round")
        else:
            cv.line(jx, y1 + 24, jx, y2 - 24, stroke=LINE, stroke_width=2.5,
                    opacity=.5, cap="round")
        cv.circle(jx, y1 + 24, 3.4, fill=MUT)
        cv.circle(jx, y2 - 24, 3.4, fill=MUT)

    # 5. merged wells
    for pr in merged:
        a = pos.get(pr[0])
        b = pos.get(pr[1])
        if not a or not b:
            continue
        x = min(a[0], b[0]) - 20
        y = min(a[1], b[1]) - 19
        w = abs(a[0] - b[0]) + 40
        h = abs(a[1] - b[1]) + 38
        cv.rect(x, y, w, h, rx=11, fill=PURPLE, fill_opacity=.10,
                stroke=PURPLE, stroke_width=1.8, dash=[5, 3])
        cv.text(x + w / 2, y - 5, "well", 9.5, PURPLE, anchor="middle",
                weight="600")

    # 6. ions
    for iid, (lab, typ) in ions.items():
        p = pos.get(iid)
        if not p:
            continue
        cv.push()
        cv.translate(p[0], p[1])
        if typ == "data":
            cv.circle(0, 0, 15, fill=PANEL, stroke=LINE, stroke_width=1.4)
            cv.text(0, 4, lab, 10, INK, anchor="middle")
        else:
            col = {"X": C_X, "Z": C_Z, "comm": TEAL, "spare": LINE}[typ]
            dk = {"X": C_XD, "Z": C_ZD, "comm": COMM_DK, "spare": MUT}[typ]
            on = iid in hi
            cv.rect(-14, -13, 28, 26, rx=5, fill=col,
                    stroke=(AMBER if on else dk),
                    stroke_width=(2.8 if on else 1.0))
            tcol = MUT if typ == "spare" else "#ffffff"
            cv.text(0, 4, lab, 9, tcol, anchor="middle", weight="600")
        cv.pop()


def build_data():
    FR, ions, home = V.build(merge=True, rounds=1, d=3)
    celly = [V.CY[r] for r in range(V.D)]
    junctions = [{"c": c, "b": b, "x": V.JX(c), "y1": V.CY[b], "y2": V.CY[b + 1]}
                 for c in range(V.D) for b in range(V.D - 1)]
    npark = sum(1 for s in V.STABS if V.is_right_boundary(s))
    wells = ([home[i] for i in ions if ions[i][1] in ("data", "X", "Z")]
             + [[V.X(V.D - 0.5) + 48 * (k + 1), V.CY[V.D - 1]]
                for k in range(npark)])
    xlo = V.X(-0.5) - 45
    # +40 over the visualizer's XIF+70 so the "interface (cavities)" label and
    # the teal cavity arcs sit fully inside the panel instead of clipping at the
    # right edge (the live SVG clips them; a static thesis figure should not).
    xhi = V.XIF + 110
    xif = V.XIF
    DATA = {"ions": ions, "frames": FR, "celly": celly, "junctions": junctions,
            "wells": wells, "xlo": xlo, "xhi": xhi, "xif": xif}
    return DATA, FR


SEL = [5, 10, 11, 18]
CAPS = ["(a) In-row gate", "(b) Cross-row lift", "(c) Cross-row gate",
        "(d) Seam"]


def compose(DATA, FR, out_pdf, out_png=None, dpi=140):
    W = DATA["xhi"]
    H = DATA["celly"][-1] + 70
    PANEL_W = W
    PANEL_H = H
    CAP_H = 26.0
    PAD_X = 18.0
    GAP_Y = 14.0
    PAD_TOP = 10.0
    PAD_BOT = 12.0

    page_w = PANEL_W + 2 * PAD_X
    page_h = (PAD_TOP + 4 * (CAP_H + PANEL_H) + 3 * GAP_Y + PAD_BOT)

    target_w_pt = 7.0 * 72.0
    s_pt = target_w_pt / page_w
    surf_pdf = cairo.PDFSurface(out_pdf, page_w * s_pt, page_h * s_pt)
    ctx = cairo.Context(surf_pdf)
    ctx.scale(s_pt, s_pt)
    _paint(ctx, DATA, FR, page_w, page_h, PANEL_W, PANEL_H,
           CAP_H, PAD_X, GAP_Y, PAD_TOP)
    surf_pdf.finish()

    if out_png:
        s_px = dpi / 72.0 * s_pt
        pw = int(round(page_w * s_px))
        ph = int(round(page_h * s_px))
        surf_png = cairo.ImageSurface(cairo.FORMAT_ARGB32, pw, ph)
        cpx = cairo.Context(surf_png)
        cpx.scale(s_px, s_px)
        _paint(cpx, DATA, FR, page_w, page_h, PANEL_W, PANEL_H,
               CAP_H, PAD_X, GAP_Y, PAD_TOP)
        surf_png.write_to_png(out_png)


def _paint(ctx, DATA, FR, page_w, page_h, PANEL_W, PANEL_H, CAP_H, PAD_X,
           GAP_Y, PAD_TOP):
    cv = Canvas(ctx)
    cv.rect(0, 0, page_w, page_h, fill="#ffffff")
    y = PAD_TOP
    for k, cap in zip(SEL, CAPS):
        cv.text(PAD_X + 2, y + 17, cap, 15, INK, weight="bold")
        y += CAP_H
        cv.rect(PAD_X, y, PANEL_W, PANEL_H, rx=10, fill=BG, stroke=LINE,
                stroke_width=1.0)
        cv.push()
        cv.translate(PAD_X, y)
        draw_frame(cv, DATA, FR[k])
        cv.pop()
        y += PANEL_H + GAP_Y


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    DATA, FR = build_data()
    for k, cap in zip(SEL, CAPS):
        print("frame %2d  %-22s  %s" % (k, cap, FR[k]["cap"][:60]))
    out_pdf = os.path.join(here, "ch4_round.pdf")
    out_png = os.path.join(here, "ch4_round_preview.png")
    compose(DATA, FR, out_pdf, out_png, dpi=140)
    print("wrote", out_pdf)
    print("wrote", out_png)


if __name__ == "__main__":
    main()

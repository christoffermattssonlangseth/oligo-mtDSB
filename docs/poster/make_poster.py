#!/usr/bin/env python
"""Generate an editable A0 poster (docs/poster/poster.svg) for the OL responder finding.

Karolinska Institutet style, after KaroSpace/poster/karospace-poster.svg: plum gradient header with a
coral strip, light plum background, white cards with numbered plum title bands, big-number stats row.
Data glyphs keep the figure colours (blue = control, orange = responder, open orange = non-responder).

Plain SVG: every text block is a <text> element, every panel a <g id=...> group, and figures are
*linked* (relative paths into results/), so the file stays small and picks up regenerated figures.
For print, vector PDFs of every figure sit next to the PNGs in results/*/figures/.

    python docs/poster/make_poster.py
    rsvg-convert -w 2400 docs/poster/poster.svg -o docs/poster/poster_preview.png
"""

from __future__ import annotations

import random
from pathlib import Path
from xml.sax.saxutils import escape

HERE = Path(__file__).resolve().parent
W, H = 841, 1189  # A0 portrait, mm
M = 25            # outer margin
GAP = 12

# Karolinska palette (as in the KaroSpace poster)
PLUM, PLUM_DARK, CORAL = "#870052", "#4f0433", "#ff876f"
BG, CARD, BORDER = "#f6f0f4", "#ffffff", "#e3d6de"
INK, MUTED, PINK_LIGHT, PINK_MID = "#1a1a1a", "#6a5667", "#f6dcec", "#d7aec6"
# data colours (match the figures)
BLUE, ORANGE, GREYDOT = "#2a78d6", "#eb6834", "#a0889a"
FONT = "'Helvetica Neue', Helvetica, Arial, sans-serif"

out: list[str] = []


def add(s: str) -> None:
    out.append(s)


def wrap(text: str, width_mm: float, size: float) -> list[str]:
    """Greedy wrap using an average glyph width of 0.5 em (Helvetica)."""
    max_chars = max(8, int(width_mm / (0.5 * size)))
    lines, cur = [], ""
    for word in text.split():
        if len(cur) + len(word) + (1 if cur else 0) > max_chars:
            lines.append(cur)
            cur = word
        else:
            cur = f"{cur} {word}" if cur else word
    if cur:
        lines.append(cur)
    return lines


def text(x, y, s, size=8.5, color=INK, weight="normal", anchor="start", style="normal") -> None:
    add(f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" fill="{color}" '
        f'font-weight="{weight}" font-style="{style}" text-anchor="{anchor}">{escape(s)}</text>')


def para(x, y, s, width, size=8.5, color=INK, lh=1.3, weight="normal") -> float:
    lines = wrap(s, width, size)
    add(f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" fill="{color}" font-weight="{weight}">')
    for i, ln in enumerate(lines):
        add(f'<tspan x="{x:.1f}" dy="{0 if i == 0 else size * lh:.2f}">{escape(ln)}</tspan>')
    add("</text>")
    return y + size * lh * len(lines)


def bullets(x, y, items, width, size=8.5, color=INK, lh=1.3, gap=2.5, mark=PLUM) -> float:
    for it in items:
        text(x, y, "✓", size=size, color=mark, weight="bold")
        y = para(x + size * 0.9, y, it, width - size * 0.9, size, color, lh) + gap
    return y


def panel(gid, x, y, w, h, title, band=PLUM, title_size=11):
    """White card with a plum title band; returns content origin and width."""
    add(f'<g id="{gid}">')
    add(f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="3" fill="{CARD}" stroke="{BORDER}" stroke-width="0.8"/>')
    add(f'<path d="M{x + 3:.1f},{y:.1f} H{x + w - 3:.1f} Q{x + w:.1f},{y:.1f} {x + w:.1f},{y + 3:.1f} V{y + 17:.1f} '
        f'H{x:.1f} V{y + 3:.1f} Q{x:.1f},{y:.1f} {x + 3:.1f},{y:.1f} Z" fill="{band}"/>')
    text(x + 8, y + 12, title, size=title_size, color="#ffffff", weight="bold")
    return x + 8, y + 28, w - 16


def end():
    add("</g>")


def image(x, y, w, h, href, gid=None):
    idattr = f' id="{gid}"' if gid else ""
    add(f'<image{idattr} x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" href="{href}" '
        f'xlink:href="{href}" preserveAspectRatio="xMidYMid meet"/>')


def mouse(cx, cy, cond, resp=None, r=4.2):
    if cond == "control":
        add(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r}" fill="{BLUE}"/>')
    elif resp:
        add(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r}" fill="{ORANGE}"/>')
    else:
        add(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r - 0.8}" fill="#ffffff" stroke="{ORANGE}" stroke-width="1.6"/>')


def arrow(x1, y1, x2, y2, color=MUTED, w=1.6):
    add(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{color}" stroke-width="{w}" '
        f'marker-end="url(#arrow)"/>')


# figures are copied next to the poster so the folder is self-contained (and renderers that refuse
# links outside the SVG's folder still show them); rerun this script after regenerating figures
import shutil
REPO = HERE.parents[1]
(HERE / "figures").mkdir(exist_ok=True)
for src in ["results/sample_classifier/figures/F2_three_readouts", "results/sample_classifier/figures/F9_ol_isr_maps",
            "results/ol_gene_drivers/figures/G3_driver_genes_by_section"]:
    for ext in (".png", ".pdf"):
        shutil.copy2(REPO / (src + ext), HERE / "figures" / (Path(src).name + ext))
FIG = GEN = "figures/"

# ----------------------------------------------------------------------------------------------
add('<?xml version="1.0" encoding="UTF-8"?>')
add(f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
    f'width="{W}mm" height="{H}mm" viewBox="0 0 {W} {H}">')
add("<title>oligo-mtDSB poster</title>")
add(f'''<defs>
  <linearGradient id="header-gradient" x1="0" y1="0" x2="0.35" y2="1">
    <stop offset="0" stop-color="{PLUM_DARK}"/><stop offset="1" stop-color="{PLUM}"/></linearGradient>
  <linearGradient id="footer-gradient" x1="0" y1="0" x2="1" y2="0">
    <stop offset="0" stop-color="{PLUM}"/><stop offset="1" stop-color="{PLUM_DARK}"/></linearGradient>
  <marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
    <path d="M0,0 L10,5 L0,10 z" fill="{MUTED}"/></marker>
</defs>''')
add(f'<rect id="background" width="{W}" height="{H}" fill="{BG}"/>')

# ---- header ----------------------------------------------------------------------------------
HH = 104
add('<g id="header">')
add(f'<rect x="0" y="0" width="{W}" height="{HH}" fill="url(#header-gradient)"/>')
add(f'<rect x="0" y="{HH - 3.5}" width="{W}" height="3.5" fill="{CORAL}"/>')
text(M, 36, "Mitochondrial DNA breaks trigger an all-or-none", size=24, color="#ffffff", weight="bold")
text(M, 63, "stress response in oligodendrocytes", size=24, color="#ffffff", weight="bold")
text(M, 79, "Leave-one-mouse-out classifiers on Xenium 5k spatial transcriptomics", size=10.5, color=PINK_LIGHT)
text(M, 93, "Christoffer M. Langseth, [co-authors]  ·  Karolinska Institutet, [affiliations]", size=7.5, color=PINK_MID)
image(W - M - 82, 20, 82, 82 * 556 / 1356, "ki-logo-white.png", gid="ki_logo")
add("</g>")

# ---- stats row -------------------------------------------------------------------------------
yS = HH + GAP
hS = 40
add('<g id="stats_row">')
stats = [("3 / 6", "mtDSB mice respond", PLUM), ("43", "driver genes (34 up, 9 down)", CORAL),
         ("10 / 50", "top genes replicate across cohorts", PLUM_DARK), ("22 / 22", "agree with earlier Bayesian model", PLUM)]
sw = (W - 2 * M) / 4
for i, (num, lab, col) in enumerate(stats):
    sx = M + i * sw + sw / 2
    text(sx, yS + 22, num, size=22, color=col, weight="bold", anchor="middle")
    text(sx, yS + 33, lab, size=7.5, color=MUTED, anchor="middle")
add("</g>")

# ---- row A: background | design | earlier analyses -------------------------------------------
yA, hA = yS + hS + GAP, 142
cw = (W - 2 * M - 2 * GAP) / 3
xs = [M, M + cw + GAP, M + 2 * (cw + GAP)]

cx, cy, cwid = panel("background_panel", xs[0], yA, cw, hA, "1 · Background")
cy = para(cx, cy, "Mitochondrial dysfunction in oligodendrocytes (OLs) is implicated in white-matter disease. "
          "In PlpTA:mtPst1 mice, removing doxycycline induces double-strand breaks in OL mitochondrial DNA (mtDSB).",
          cwid, size=7.8)
cy += 3
cy = para(cx, cy, "This line shows demyelination of corpus callosum and cerebellum 4-6 months after induction. "
          "What happens in the tissue before that, and in which cells and regions?", cwid, size=7.8)
cy += 3
para(cx, cy, "Can the tissue itself tell which group a mouse belongs to, and which genes does it use?",
     cwid, size=7.8, weight="bold", color=PLUM)
end()

cx, cy, cwid = panel("design_panel", xs[1], yA, cw, hA, "2 · Study design")
ty = cy + 6
add(f'<line x1="{cx}" y1="{ty}" x2="{cx + cwid}" y2="{ty}" stroke="{MUTED}" stroke-width="1"/>')
for frac, lab, sub, anc in [(0.0, "DOX removal", "P21 or P60", "start"), (0.5, "+42 days", "", "middle"),
                            (1.0, "Xenium 5k", "1 section / mouse", "end")]:
    px = cx + frac * cwid
    add(f'<circle cx="{px:.1f}" cy="{ty}" r="2.4" fill="{PLUM}"/>')
    text(px, ty - 5, lab, size=7.2, color=INK, anchor=anc, weight="bold")
    if sub:
        text(px, ty + 8.5, sub, size=6.5, color=MUTED, anchor=anc)
gy = ty + 30
ids = {"21": [("RB4620", "control", None), ("RB4658", "control", None), ("RB4676", "control", None),
              ("RB4627", "mtDSB", True), ("RB4630", "mtDSB", False), ("RB4653", "mtDSB", False)],
       "60": [("RB4282", "control", None), ("RB4350", "control", None), ("RB4498", "control", None),
              ("RB4401", "mtDSB", True), ("RB4405", "mtDSB", True), ("RB4403", "mtDSB", False)]}
colx = [cx + 50 + i * 30 for i in range(6)]
text(colx[1], gy - 9, "control", size=7.2, color=BLUE, anchor="middle", weight="bold")
text(colx[4], gy - 9, "mtDSB", size=7.2, color=ORANGE, anchor="middle", weight="bold")
for row, age in enumerate(["21", "60"]):
    yy = gy + row * 25
    text(cx, yy + 2.5, f"P{age} cohort", size=7.2, color=MUTED)
    for (sid, cond, resp), xx in zip(ids[age], colx):
        mouse(xx, yy, cond, resp)
        text(xx, yy + 11, sid[2:], size=5.6, color=MUTED, anchor="middle")
ly = gy + 58
mouse(cx + 3, ly, "mtDSB", True, r=3); text(cx + 9, ly + 2.4, "responder", size=6.5, color=MUTED)
mouse(cx + 50, ly, "mtDSB", False, r=3); text(cx + 56, ly + 2.4, "non-responder", size=6.5, color=MUTED)
mouse(cx + 108, ly, "control", r=3); text(cx + 114, ly + 2.4, "control", size=6.5, color=MUTED)
para(cx, ly + 13, "980 474 cells · 5 101 genes · 14 cell classes · 11 brain regions", cwid, size=6.8, color=MUTED)
end()

cx, cy, cwid = panel("earlier_panel", xs[2], yA, cw, hA, "3 · Earlier analyses, explained")
cy = para(cx, cy, "Group comparisons (DESeq2, hierarchical Bayes) gave weak, age-dependent effects. "
          "Averaging 3 responders with 3 non-responders explains both:", cwid, size=7.6)
cy += 3
bullets(cx, cy, [
    "Bayesian case-control model: all 22 driver genes it tested agree in direction (p < 1e-6), but shifts were small (0.01-0.05).",
    "'P21 modest / P60 strong' = 1 of 3 vs 2 of 3 responders; the P21 responder is as strong as any.",
    "Top 'mtDSB genes' Sgk1, Fkbp5 (stress hormones) varied between sections, not with treatment.",
], cwid, size=7.3, gap=1.8)
end()

# ---- methods ---------------------------------------------------------------------------------
yB, hB = yA + hA + GAP, 196
cx, cy, cwid = panel("methods_panel", M, yB, W - 2 * M, hB, "4 · Methods: let the tissue predict the mouse",
                     band=PLUM_DARK, title_size=12)
n = 5
bw = (cwid - (n - 1) * 14) / n
by, bh = cy, 124
steps = [
    ("Spatial data", "Xenium 5k, 12 mice. Counts summed per cell type x region x mouse (pseudobulk) or kept per cell."),
    ("General labels", "Disease-defined cluster (DA-OL) folded into oligodendrocytes, so the model, not the clustering, must find the effect."),
    ("Leave-one-mouse-out", "Train on 11 mice, predict the 12th; repeat for every mouse. Ridge logistic regression; permutation tests at the mouse level."),
    ("Three readouts", "Prespecified ISR gene score; classifier score per OL; DA-OL cluster share. Agreement between independent inputs."),
    ("Driver genes", "Correlate each gene with the OL score within each section (removes batch, plane, sacrifice time). Replicate across cohorts."),
]
for i, (head, body) in enumerate(steps):
    bx = cx + i * (bw + 14)
    add(f'<g id="method_step_{i + 1}">')
    add(f'<rect x="{bx:.1f}" y="{by:.1f}" width="{bw:.1f}" height="{bh}" rx="3" fill="{BG}"/>')
    add(f'<circle cx="{bx + 9:.1f}" cy="{by + 9:.1f}" r="5.5" fill="{CORAL if i in (2, 4) else PLUM}"/>')
    text(bx + 9, by + 11.6, str(i + 1), size=7.5, color="#ffffff", weight="bold", anchor="middle")
    text(bx + 18, by + 11.5, head, size=8.3, color=PLUM_DARK, weight="bold")
    para(bx + 5, by + 25, body, bw - 10, size=7.0, color=INK)
    iy = by + bh - 32
    if i == 0:
        for k in range(18):
            add(f'<circle cx="{bx + 12 + (k % 9) * 13:.1f}" cy="{iy + (k // 9) * 12:.1f}" r="3.2" fill="{BORDER if k % 4 else PLUM}"/>')
    elif i == 1:
        add(f'<circle cx="{bx + 22:.1f}" cy="{iy + 5:.1f}" r="8" fill="{PINK_LIGHT}"/>')
        add(f'<circle cx="{bx + 40:.1f}" cy="{iy + 5:.1f}" r="6" fill="{ORANGE}" opacity="0.8"/>')
        arrow(bx + 52, iy + 5, bx + 74, iy + 5)
        add(f'<circle cx="{bx + 92:.1f}" cy="{iy + 5:.1f}" r="11" fill="{PINK_LIGHT}"/>')
        text(bx + 92, iy + 25, "Oligodendrocytes", size=6.0, color=MUTED, anchor="middle")
    elif i == 2:
        for k in range(12):
            xx, yy = bx + 12 + (k % 6) * 13, iy + (k // 6) * 13
            if k == 9:
                add(f'<circle cx="{xx:.1f}" cy="{yy:.1f}" r="4" fill="#ffffff" stroke="{CORAL}" stroke-width="1.8"/>')
            else:
                add(f'<circle cx="{xx:.1f}" cy="{yy:.1f}" r="3.6" fill="{GREYDOT}"/>')
        text(bx + 12, iy + 26, "11 train", size=6.0, color=MUTED)
        text(bx + 52, iy + 26, "1 predict", size=6.0, color=CORAL, weight="bold")
    elif i == 3:
        for k, (lab, col) in enumerate([("ISR genes (prior)", PLUM_DARK), ("classifier (label-free)", PLUM), ("cluster (descriptive)", CORAL)]):
            add(f'<rect x="{bx + 8:.1f}" y="{iy - 4 + k * 10:.1f}" width="6" height="6" rx="1" fill="{col}"/>')
            text(bx + 18, iy + 1.5 + k * 10, lab, size=6.5, color=MUTED)
    elif i == 4:
        random.seed(1)
        x0, y0 = bx + 10, iy + 22
        add(f'<line x1="{x0}" y1="{y0}" x2="{x0 + 60}" y2="{y0}" stroke="{BORDER}"/>')
        add(f'<line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y0 - 26}" stroke="{BORDER}"/>')
        for _ in range(28):
            t = random.random()
            add(f'<circle cx="{x0 + 4 + t * 54:.1f}" cy="{y0 - 3 - t * 20 + random.gauss(0, 3):.1f}" r="1.3" fill="{PLUM}"/>')
        text(x0 + 66, y0 - 12, "gene vs", size=6.0, color=MUTED)
        text(x0 + 66, y0 - 5, "OL score", size=6.0, color=MUTED)
    add("</g>")
    if i < n - 1:
        arrow(bx + bw + 1, by + bh / 2, bx + bw + 13, by + bh / 2)
sy = by + bh + 14
text(cx, sy, "Why this is robust", size=9.5, color=PLUM, weight="bold")
rob = ["The mouse is the unit: no cell is ever scored by a model that saw its mouse.",
       "Confounds cannot create driver genes: they must vary between cells of the same section.",
       "Independent of clustering: prior gene set, label-free classifier and cluster agree.",
       "Replication: models trained on disjoint mice (P21-only vs P60-only) share top genes."]
half = (cwid - 20) / 2
bullets(cx, sy + 12, rob[:2], half, size=7.5, gap=1.5)
bullets(cx + half + 20, sy + 12, rob[2:], half, size=7.5, gap=1.5)
end()

# ---- result 1: responders --------------------------------------------------------------------
yC, hC = yB + hB + GAP, 212
cx, cy, cwid = panel("result1_panel", M, yC, W - 2 * M, hC,
                     "5 · Only 3 of 6 mtDSB mice respond, and they respond strongly", title_size=12)
iw = 540
ih = iw / 3.233
image(cx, cy - 2, iw, ih, FIG + "F2_three_readouts.png", gid="fig_F2")
tx = cx + iw + 14
twid = cwid - iw - 14
add('<g id="result1_text">')
text(tx, cy + 16, "3 / 6", size=28, color=CORAL, weight="bold")
text(tx + 64, cy + 8, "mtDSB mice with an", size=7.6, color=MUTED)
text(tx + 64, cy + 17, "OL stress response", size=7.6, color=MUTED)
bullets(tx, cy + 34, [
    "Same 3 mice by three independent readouts.",
    "ISR genes (Atf5, Trib3, Ddit3, Sesn2): no overlap with controls; non-responders = controls.",
    "mtDSB group 10x more variable than controls (p = 1e-4) while means do not differ.",
    "Not explained by sex, genotype, sample swap, white-matter content or Xenium slide.",
], twid, size=7.3, gap=1.4)
add("</g>")
para(cx, cy + ih + 6, "Each mouse scored by models that never saw it. Filled orange = responders, open = non-responders, blue = controls.",
     iw, size=6.5, color=MUTED)
end()

# ---- row D: genes | where --------------------------------------------------------------------
yD = yC + hC + GAP
hD = 336
lw = 262
cx, cy, cwid = panel("result2_panel", M, yD, lw, hD, "6 · What responding OLs do", title_size=11)
gh = 236
image(cx + (cwid - gh * 0.851) / 2, cy - 4, gh * 0.851, gh, GEN + "G3_driver_genes_by_section.png", gid="fig_G3")
bullets(cx, cy + gh + 2, [
    "34 up: ATF4 stress response, mitochondrial & protein quality control (Hspa9, p62), p21, MHC-I and Cd40.",
    "9 myelin genes down (Mog, Cnp, Cldn11, Mal).",
    "Top genes replicate across cohorts (10/50 shared, p = 2e-11).",
], cwid, size=7.0, gap=1.1)
end()

rx = M + lw + GAP
rw = W - M - rx
cx, cy, cwid = panel("result3_panel", rx, yD, rw, hD, "7 · Where: the stress response in space (ISR score)",
                     title_size=11)
ih3 = cwid / (3139 / 1481)
image(cx, cy - 6, cwid, ih3, FIG + "F9_ol_isr_maps.png", gid="fig_F9")
ty3 = cy + ih3 + 2
half3 = (cwid - 16) / 2
bullets(cx, ty3, [
    "Prespecified ISR score in every OL (spatially smoothed, percentile of same-slide controls): 31-55% of OLs above the control 90th percentile in responders vs 7-14% elsewhere.",
    "Blue = below, red = above the control median. In controls, white-matter OLs have the lowest ISR; in responders the thalamus and fiber-tract arcs turn red.",
], half3, size=7.3, gap=1.3)
bullets(cx + half3 + 16, ty3, [
    "Thalamus and fiber tracts rank top on all three readouts; myelin genes drop most there.",
    "Microglia in the same three mice turn DAM-like (Cd74, Cst7 10-30x).",
], half3, size=7.3, gap=1.3)
end()

# ---- footer: conclusions --------------------------------------------------------------------
yE = yD + hD + GAP
hE = H - yE
add('<g id="conclusions_panel">')
add(f'<rect x="0" y="{yE - 3.5:.1f}" width="{W}" height="3.5" fill="{CORAL}"/>')
add(f'<rect x="0" y="{yE:.1f}" width="{W}" height="{hE:.1f}" fill="url(#footer-gradient)"/>')
text(M, yE + 15, "8 · Conclusions and next steps", size=11, color="#ffffff", weight="bold")
half = (W - 2 * M - 20) / 2
bullets(M, yE + 28, [
    "mtDSB causes an all-or-none OL stress response at 42 days: ISR + mitochondrial quality control, immune visibility, myelin loss.",
    "Responding OLs sit with DAM-like microglia in thalamus and white matter: a candidate prelude to demyelination.",
], half, size=7.4, color="#ffffff", mark=CORAL, gap=1.4)
bullets(M + half + 20, yE + 28, [
    "Why do 3 mice not respond? Testing induction (mtPstI qPCR, mtDNA damage) on remaining tissue.",
    "Code, figures, tables: github.com/christoffermattssonlangseth/oligo-mtDSB",
], half, size=7.4, color="#ffffff", mark=PINK_LIGHT, gap=1.4)
add("</g>")

add("</svg>")
(HERE / "poster.svg").write_text("\n".join(out))
print("wrote", HERE / "poster.svg")

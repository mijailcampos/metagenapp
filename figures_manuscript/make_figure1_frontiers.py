#!/usr/bin/env python3
"""Figure 1 of the manuscript in Frontiers format (MetagenApp v1.1 workflow, 9 steps +
downstream analyses), in two rows at the final two-column width (180 mm) with all text >= 8 pt.
Saved as RGB TIFF (LZW, 300 dpi) plus a PNG preview in output/."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Ellipse, Polygon, Circle
from matplotlib.lines import Line2D
from PIL import Image

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial"]})
HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
os.makedirs(HERE, exist_ok=True)
W = 180 / 25.4
TINT = {"#2a78d6": "#eaf2fc", "#eb6834": "#fdece3", "#1baf7a": "#e7f8f1", "#eda100": "#fdf2dc"}
INK = "#33322f"
steps = [
    ("1. Input", "Paired-end\nFASTQ reads", "#2a78d6"),
    ("2. Assembly", "VSEARCH merge\n(≤ 5 mismatches)", "#2a78d6"),
    ("3. Primer\ntrimming", "cutadapt\n341F / 805R", "#eb6834"),
    ("4. Quality\nfiltering", "400–440 bp, N = 0,\nhomopolymer ≤ 8", "#eb6834"),
    ("5. Dereplication", "Unique sequences\nwith abundances", "#1baf7a"),
    ("6. Chimera\nremoval", "UCHIME de novo\n(abundance-aware)", "#1baf7a"),
    ("7. Reference\nfiltering", "VSEARCH vs SILVA\n(identity ≥ 0.70)", "#1baf7a"),
    ("8. Feature\nclustering", "97% OTUs (UPARSE)\nsingletons mapped", "#eda100"),
    ("9. Taxonomic\nclassification", "Wang classifier\nSILVA 138.2 NR99", "#eda100"),
]
NCOL = 5
box_w, box_h = 0.92, 1.15
row_y = [1.55, 0.0]  # bottom y of row 1 and row 2
fig_w_units = NCOL
fig_h_units = row_y[0] + box_h + 0.05
fig = plt.figure(figsize=(W, W * fig_h_units / fig_w_units))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, fig_w_units)
ax.set_ylim(-0.05, fig_h_units)
ax.axis("off")
# ---- per-step icon glyphs (drawn centered at (cx, cy), ~0.3 x 0.3 footprint) ----

def icon_input(cx, cy, color):
    for off, alpha in ((0.045, 1.0), (-0.045, 0.55)):
        ax.add_patch(FancyBboxPatch((cx - 0.09, cy - 0.11 + off), 0.15, 0.16,
                                     boxstyle="round,pad=0.006,rounding_size=0.015",
                                     linewidth=1.1, edgecolor=color, facecolor="white",
                                     alpha=alpha, zorder=3))
    for ly in (-0.03, 0.01, 0.05):
        ax.add_line(Line2D([cx - 0.06, cx + 0.03], [cy + 0.045 + ly, cy + 0.045 + ly],
                            color=color, linewidth=0.8, alpha=0.8, zorder=4))

def icon_assembly(cx, cy, color):
    ax.add_patch(FancyArrowPatch((cx - 0.14, cy + 0.06), (cx - 0.01, cy + 0.06),
                                  arrowstyle="-|>", mutation_scale=8, linewidth=1.3, color=color))
    ax.add_patch(FancyArrowPatch((cx + 0.14, cy - 0.06), (cx + 0.01, cy - 0.06),
                                  arrowstyle="-|>", mutation_scale=8, linewidth=1.3, color=color))
    ax.add_patch(Rectangle((cx - 0.09, cy - 0.10), 0.18, 0.20, fill=False,
                            edgecolor=color, linewidth=1.0, linestyle=(0, (1.5, 1.2))))

def icon_primer(cx, cy, color):
    for sx in (-1, 1):
        ax.add_line(Line2D([cx, cx + sx * 0.11], [cy, cy + 0.11], color=color, linewidth=1.3))
        ax.add_line(Line2D([cx, cx + sx * 0.11], [cy, cy - 0.11], color=color, linewidth=1.3))
        ax.add_patch(Circle((cx + sx * 0.11, cy + 0.11), 0.022, facecolor="white",
                             edgecolor=color, linewidth=1.0))
        ax.add_patch(Circle((cx + sx * 0.11, cy - 0.11), 0.022, facecolor="white",
                             edgecolor=color, linewidth=1.0))
    ax.add_patch(Circle((cx, cy), 0.018, facecolor=color, edgecolor="none"))

def icon_quality(cx, cy, color):
    ax.add_patch(Polygon([(cx - 0.13, cy + 0.11), (cx + 0.13, cy + 0.11),
                           (cx + 0.03, cy - 0.02), (cx + 0.03, cy - 0.11),
                           (cx - 0.03, cy - 0.11), (cx - 0.03, cy - 0.02)],
                          closed=True, facecolor="white", edgecolor=color, linewidth=1.2))

def icon_derep(cx, cy, color):
    for off, alpha in ((0.05, 0.4), (0.025, 0.65), (0.0, 1.0)):
        ax.add_patch(Rectangle((cx - 0.09 + off, cy - 0.09 + off), 0.14, 0.14,
                                facecolor="white", edgecolor=color, linewidth=1.0, alpha=alpha))

def icon_chimera(cx, cy, color):
    ax.add_line(Line2D([cx - 0.13, cx - 0.03, cx + 0.02, cx + 0.13],
                        [cy - 0.04, cy + 0.08, cy - 0.06, cy + 0.06],
                        color=color, linewidth=1.3, alpha=0.55))
    for sx, sy in ((-1, -1), (-1, 1), (1, -1), (1, 1)):
        ax.add_line(Line2D([cx + sx * 0.03, cx + sx * 0.10], [cy + sy * 0.03, cy + sy * 0.10],
                            color="#c0392b", linewidth=1.6))

def icon_reffilter(cx, cy, color):
    # query bar aligned to a reference bar, with a check mark
    ax.add_patch(Rectangle((cx - 0.13, cy + 0.04), 0.26, 0.045, facecolor=color, edgecolor="none", alpha=0.45))
    ax.add_patch(Rectangle((cx - 0.08, cy - 0.03), 0.16, 0.045, facecolor="white", edgecolor=color, linewidth=1.1))
    ax.add_line(Line2D([cx - 0.04, cx - 0.005, cx + 0.07], [cy - 0.09, cy - 0.125, cy - 0.06],
                        color=color, linewidth=1.6))

def icon_cluster(cx, cy, color):
    ax.add_patch(Ellipse((cx, cy), 0.30, 0.24, fill=False, edgecolor=color,
                          linewidth=1.0, linestyle=(0, (1.5, 1.2))))
    for ox, oy in ((-0.06, 0.03), (0.05, 0.04), (0.0, -0.05), (0.08, -0.03)):
        ax.add_patch(Circle((cx + ox, cy + oy), 0.028, facecolor=color, edgecolor="none"))

def icon_classify(cx, cy, color):
    ax.add_patch(Circle((cx - 0.02, cy + 0.03), 0.08, fill=False, edgecolor=color, linewidth=1.4))
    ax.add_line(Line2D([cx + 0.04, cx + 0.13], [cy - 0.04, cy - 0.13], color=color, linewidth=1.8))

step_icons = [icon_input, icon_assembly, icon_primer, icon_quality,
              icon_derep, icon_chimera, icon_reffilter, icon_cluster, icon_classify]

def card(x, y0, color, fill):
    ax.add_patch(FancyBboxPatch((x - box_w / 2, y0), box_w, box_h, boxstyle="round,pad=0.02,rounding_size=0.06",
                                linewidth=1.3, edgecolor=color, facecolor=fill))

def arrow(p0, p1):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=10, linewidth=1.2, color="#898781"))

pos = [(i % NCOL + 0.5, row_y[i // NCOL]) for i in range(len(steps))]
for i, ((title, sub, color), draw) in enumerate(zip(steps, step_icons)):
    x, y0 = pos[i]
    card(x, y0, color, TINT[color])
    draw(x, y0 + box_h - 0.27, color)
    ax.text(x, y0 + box_h - 0.5, title, ha="center", va="top", fontsize=8.5, fontweight="bold", color=INK,
            linespacing=1.1)
    ax.text(x, y0 + 0.25, sub, ha="center", va="center", fontsize=8, color="#52514e", linespacing=1.15)

# arrows within rows
for i in range(len(steps) - 1):
    (xa, ya), (xb, yb) = pos[i], pos[i + 1]
    if ya == yb:
        arrow((xa + box_w / 2, ya + box_h / 2), (xb - box_w / 2, yb + box_h / 2))
# row 1 -> row 2: elbow from the bottom of step 5 to the top of step 6
xa, ya = pos[NCOL - 1]
xb, yb = pos[NCOL]
mid = (ya + yb + box_h) / 2
ax.add_line(Line2D([xa, xa, xb, xb], [ya - 0.02, mid, mid, yb + box_h + 0.12], color="#898781", linewidth=1.2))
arrow((xb, yb + box_h + 0.14), (xb, yb + box_h + 0.02))

# 10. Downstream analyses (last slot of row 2)
DOWN, DOWN_TINT = "#4a3aa7", "#efecfa"
xd, yd = NCOL - 0.5, row_y[1]
card(xd, yd, DOWN, DOWN_TINT)
ax.text(xd, yd + box_h - 0.1, "10. Downstream\nanalyses", ha="center", va="top", fontsize=8.5, fontweight="bold",
        color="#2c2358", linespacing=1.1)
for k, lab in enumerate(["Taxonomic profiles", "Alpha diversity", "Beta diversity", "Exportable results"]):
    ax.text(xd - box_w / 2 + 0.05, yd + 0.66 - k * 0.17, "• " + lab, ha="left", va="center", fontsize=8,
            color="#2c2358")
x9, y9 = pos[-1]
arrow((x9 + box_w / 2, y9 + box_h / 2), (xd - box_w / 2, yd + box_h / 2))

png = os.path.join(HERE, "Figure1_pipeline_diagram.png")
fig.savefig(png, dpi=300)
im = Image.open(png).convert("RGB")
im.save(os.path.join(HERE, "Figure1_pipeline_diagram.tif"), compression="tiff_lzw", dpi=(300, 300))
sizes = [t.get_fontsize() for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip()]
print(f"Figure1: {im.size[0]}x{im.size[1]} px = {im.size[0]/300*25.4:.0f} x {im.size[1]/300*25.4:.0f} mm, smallest text {min(sizes)} pt")

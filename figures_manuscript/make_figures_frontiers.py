#!/usr/bin/env python3
"""Figures 2-6 of the manuscript in Frontiers format: final two-column width (180 mm),
every text element >= 8 pt, panel labels (A)/(B) outside the plots, saved as RGB TIFF
(LZW) at 300 dpi plus a PNG preview in output/. Data are read from results_v1.1/ and
data/; nothing is hard-coded."""

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
TABLES = os.path.join(HERE, "..", "results_v1.1", "tables")
RAREFACTION = os.path.join(HERE, "..", "results_v1.1", "rarefaction")
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "output")
os.makedirs(OUT, exist_ok=True)
W = 180 / 25.4  # final two-column width, inches
W1 = 85 / 25.4  # final one-column width, inches
MIN_PT = 8

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 8.5,
    "axes.titlesize": 9,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.edgecolor": "#c3c2b7",
    "axes.labelcolor": "#0b0b0b",
    "text.color": "#0b0b0b",
    "xtick.color": "#52514e",
    "ytick.color": "#52514e",
    "axes.grid": True,
    "grid.color": "#e1e0d9",
    "grid.linewidth": 0.6,
    "figure.facecolor": "#ffffff",
    "axes.facecolor": "#ffffff",
    "savefig.facecolor": "#ffffff",
})

PIPE_COLOR = {"MetagenApp": "#2a78d6", "mothur": "#eb6834", "QIIME 2": "#1baf7a"}
PIPES = ["MetagenApp", "mothur", "QIIME 2"]


def save(fig, name):
    """Save at exactly W wide; check the smallest font; write PNG preview and RGB TIFF."""
    sizes = [t.get_fontsize() for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip()]
    assert min(sizes) >= MIN_PT, (name, min(sizes))
    png = os.path.join(OUT, f"{name}.png")
    fig.savefig(png, dpi=300)
    plt.close(fig)
    im = Image.open(png).convert("RGB")
    im.save(os.path.join(OUT, f"{name}.tif"), compression="tiff_lzw", dpi=(300, 300))
    print(f"{name}: {im.size[0]}x{im.size[1]} px = {im.size[0] / 300 * 25.4:.0f} x "
          f"{im.size[1] / 300 * 25.4:.0f} mm, smallest text {min(sizes)} pt")


# ----------------------------------------------------------------------
# Figure 2 - Phylum-level taxonomic composition
# ----------------------------------------------------------------------
tax = pd.read_csv(f"{TABLES}/taxonomic_composition_phylum.csv")
tax_key = {"MetagenApp": "MetagenApp", "mothur": "mothur", "QIIME 2": "QIIME2"}
name_map = {
    "Bacillota": "Firmicutes",
    "Pseudomonadota": "Proteobacteria",
    "Bacteroidota": "Bacteroidetes",
    "Actinomycetota": "Actinobacteria",
    "Fusobacteriota": "Fusobacteria",
}
main_colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]

fig, ax = plt.subplots(figsize=(W1, 3.9), layout="constrained")
x = np.arange(len(PIPES))
bottoms = np.zeros(len(PIPES))
for (src, phy), col in zip(name_map.items(), main_colors):
    vals = []
    for p in PIPES:
        row = tax[(tax.pipeline == tax_key[p]) & (tax.phylum == src)]
        vals.append(float(row.mean_rel_abundance.values[0]) * 100 if len(row) else 0)
    ax.bar(x, vals, bottom=bottoms, width=0.55, color=col, label=phy, edgecolor="#ffffff", linewidth=1.2)
    bottoms += vals
other_keys = ["Patescibacteria", "Campylobacterota", "Cyanobacteriota", "Other"]
other = [tax[(tax.pipeline == tax_key[p]) & (tax.phylum.isin(other_keys))].mean_rel_abundance.sum() * 100 for p in PIPES]
uncl = [tax[(tax.pipeline == tax_key[p]) & (tax.phylum == "Unclassified")].mean_rel_abundance.sum() * 100 for p in PIPES]
ax.bar(x, other, bottom=bottoms, width=0.55, color="#d8d7d0", label="Other (<1% each)", edgecolor="#ffffff", linewidth=1.2)
bottoms += other
ax.bar(x, uncl, bottom=bottoms, width=0.55, color="#898781", label="Unclassified", edgecolor="#ffffff", linewidth=1.2)
ax.set_xticks(x)
ax.set_ylabel("Relative abundance (%)")
ax.set_ylim(0, 100)
ax.set_axisbelow(True)
ax.spines[["top", "right"]].set_visible(False)
ax.set_xticklabels(PIPES, fontweight="bold", rotation=20, ha="right")
fig.legend(loc="outside lower center", ncol=2, frameon=False, handlelength=1.0, columnspacing=1.0)
save(fig, "Figure2_taxonomic_composition")

# ----------------------------------------------------------------------
# Figure 3 - Alpha diversity boxplots (n = 36, 4,963 reads, 1,000 rarefactions)
# ----------------------------------------------------------------------
fname = {"MetagenApp": "metagenapp", "mothur": "mothur", "QIIME 2": "qiime2"}
alpha = {p: pd.read_csv(f"{RAREFACTION}/alpha_{fname[p]}_rarefied_all36.csv") for p in PIPES}
metrics = [("chao1", "Chao1"), ("shannon", "Shannon (H)"), ("invsimpson", "Inv. Simpson")]
fig, axes = plt.subplots(1, 3, figsize=(W, 2.7), layout="constrained")
rng = np.random.default_rng(19760620)
for ax, (col, title) in zip(axes, metrics):
    data = [alpha[p][col].values for p in PIPES]
    bp = ax.boxplot(data, positions=range(3), widths=0.5, patch_artist=True, showfliers=False,
                    medianprops=dict(color="#0b0b0b", linewidth=1.4),
                    whiskerprops=dict(color="#52514e", linewidth=1.0),
                    capprops=dict(color="#52514e", linewidth=1.0),
                    boxprops=dict(linewidth=1.0, edgecolor="#52514e"))
    for patch, p in zip(bp["boxes"], PIPES):
        patch.set_facecolor(PIPE_COLOR[p])
        patch.set_alpha(0.55)
    for i, p in enumerate(PIPES):
        jitter = rng.uniform(-0.12, 0.12, size=len(data[i]))
        ax.scatter(np.full(len(data[i]), i) + jitter, data[i], s=8,
                   color=PIPE_COLOR[p], alpha=0.5, edgecolor="none", zorder=3)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["MetagenApp", "mothur", "QIIME 2"], rotation=25, ha="right")
    ax.set_title(title, fontweight="bold")
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
save(fig, "Figure3_alpha_diversity")

# ----------------------------------------------------------------------
# Figure 4 - PCoA by sampling site
# ----------------------------------------------------------------------
SITE_STYLE = {
    "pharynx": {"color": "#4a3aa7", "marker": "o", "label": "Pharynx (n = 21)"},
    "nose": {"color": "#e34948", "marker": "^", "label": "Nose (n = 15)"},
}
_var = pd.read_csv(f"{TABLES}/pcoa_variance_explained.csv").replace({"pipeline": {"QIIME2": "QIIME 2"}}).set_index("pipeline")
meta = pd.read_csv(f"{DATA}/sample_site_carrier.csv")
meta["site"] = meta["site"].replace({"nasal": "nose", "nares": "nose"})
assert set(meta["site"]) == set(SITE_STYLE), set(meta["site"])


def confidence_ellipse(x, y, ax, color, n_std=2.4477):  # ~95% for bivariate normal
    cov = np.cov(x, y)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    ell = Ellipse((0, 0), width=np.sqrt(1 + pearson) * 2, height=np.sqrt(1 - pearson) * 2,
                  facecolor="none", edgecolor=color, linewidth=1.2, linestyle="--")
    tr = transforms.Affine2D().rotate_deg(45).scale(np.sqrt(cov[0, 0]) * n_std, np.sqrt(cov[1, 1]) * n_std) \
        .translate(np.mean(x), np.mean(y))
    ell.set_transform(tr + ax.transData)
    ax.add_patch(ell)


fig, axes = plt.subplots(1, 3, figsize=(W, 2.95), layout="constrained", sharex=True, sharey=True)
for ax, p in zip(axes, PIPES):
    df = pd.read_csv(f"{TABLES}/pcoa_coords_{fname[p]}.csv").merge(meta, on="sample", how="left")
    assert len(df) == 36 and df["site"].notna().all(), (p, len(df))
    if p == "QIIME 2":  # the sign of a principal coordinate is arbitrary; aligned for visual comparison
        df["PC2"] = -df["PC2"]
    for site, st in SITE_STYLE.items():
        d = df[df.site == site]
        ax.scatter(d.PC1, d.PC2, s=18, color=st["color"], marker=st["marker"],
                   alpha=0.85, edgecolor="#ffffff", linewidth=0.5, zorder=3)
        confidence_ellipse(d.PC1.values, d.PC2.values, ax, st["color"])
    ax.set_xlabel(f"PC1 ({_var.loc[p, 'PC1_pct']}%)")
    ax.set_ylabel(f"PC2 ({_var.loc[p, 'PC2_pct']}%)")
    ax.set_title(p, fontweight="bold", color=PIPE_COLOR[p])
    ax.axhline(0, color="#e1e0d9", linewidth=0.6, zorder=0)
    ax.axvline(0, color="#e1e0d9", linewidth=0.6, zorder=0)
    ax.set_xlim(-0.75, 0.9)  # same scale in all panels
    ax.set_ylim(-0.8, 0.8)
    ax.set_xticks([-0.5, 0, 0.5])
    ax.set_yticks([-0.6, -0.3, 0, 0.3, 0.6])
    ax.tick_params(labelleft=True)
    ax.spines[["top", "right"]].set_visible(False)
handles = [Line2D([0], [0], marker=st["marker"], color="none", markerfacecolor=st["color"],
                  markeredgecolor="none", markersize=6, label=st["label"]) for st in SITE_STYLE.values()]
fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False)
save(fig, "Figure4_pcoa")

# ----------------------------------------------------------------------
# Figure 5 - Genus-level taxonomic resolution (read-based | feature-based)
# ----------------------------------------------------------------------
_g = pd.read_csv(f"{TABLES}/genus_resolution_table6.tsv", sep="\t", index_col=0)
genus = {key: {p: tuple(round(float(_g.loc[f"{key2}_{c}", p]), 1)
                        for c in ["Classified", "Ambiguous/propagated", "Unclassified"]) for p in PIPES}
         for key, key2 in [("read", "reads"), ("feat", "features")]}
status_colors = {"Classified at genus level": "#2a78d6",
                 "Ambiguous / propagated": "#eda100",
                 "Unclassified at genus level": "#898781"}
fig, axes = plt.subplots(1, 2, figsize=(W, 3.3), sharey=True, layout="constrained")
for ax, key, title in zip(axes, ["read", "feat"], ["Read-based", "Feature-based"]):
    x = np.arange(3)
    c = [genus[key][p][0] for p in PIPES]
    a = [genus[key][p][1] for p in PIPES]
    u = [genus[key][p][2] for p in PIPES]
    b1 = ax.bar(x, c, width=0.55, color=status_colors["Classified at genus level"])
    b2 = ax.bar(x, a, bottom=c, width=0.55, color=status_colors["Ambiguous / propagated"])
    b3 = ax.bar(x, u, bottom=[ci + ai for ci, ai in zip(c, a)], width=0.55,
                color=status_colors["Unclassified at genus level"])
    for i in range(3):
        if c[i] > 6:
            ax.text(i, c[i] / 2, f"{c[i]:.1f}%", ha="center", va="center", color="white", fontsize=8, fontweight="bold")
        if a[i] > 4:
            ax.text(i, c[i] + a[i] / 2, f"{a[i]:.1f}%", ha="center", va="center", color="white", fontsize=8)
        if u[i] > 4:
            ax.text(i, c[i] + a[i] + u[i] / 2, f"{u[i]:.1f}%", ha="center", va="center", color="white", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(PIPES)
    ax.set_title(title, fontweight="bold")
    ax.set_ylim(0, 100)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_axisbelow(True)
axes[0].set_ylabel("Percentage (%)")
fig.legend([b1, b2, b3], list(status_colors), loc="outside lower center", ncol=3, frameon=False)
save(fig, "Figure5_genus_resolution")

# ----------------------------------------------------------------------
# Figure 6 - Mock community accuracy: (A) composition, (B) observed vs expected
# ----------------------------------------------------------------------
D = json.load(open(f"{DATA}/mock_profiles.json"))
MOCKS = [("HM-782D", "HM-782D (even)"), ("HM-783D", "HM-783D (staggered)"), ("D6306", "ZymoBIOMICS D6306")]
OTHER = "Other / unexpected"
genera = []
for m, _ in MOCKS:
    for g in sorted(D["ref"][m], key=lambda k: -D["ref"][m][k]):
        if g not in genera:
            genera.append(g)
cmap = plt.get_cmap("tab20")
GCOL = {g: cmap(i % 20) for i, g in enumerate(genera)}
GCOL[OTHER] = "#5f5e5a"

fig = plt.figure(figsize=(W, 7.6), layout="constrained")
top, bot = fig.subfigures(2, 1, height_ratios=[1.2, 1], hspace=0.04)
axs = top.subplots(1, 3)
for k, ((m, title), ax) in enumerate(zip(MOCKS, axs)):
    bars = [("Expected", D["ref"][m])] + [(p, D["obs"][f"{p}|{m}"]) for p in PIPES]
    for i, (lab, prof) in enumerate(bars):
        bottom = 0
        for g in [x for x in genera if x in D["ref"][m]] + [OTHER]:
            v = prof.get(g, 0)
            if v <= 0:
                continue
            ax.bar(i, v, bottom=bottom, width=0.72, color=GCOL[g], edgecolor="#ffffff", linewidth=0.5)
            bottom += v
    ax.set_xticks(range(4))
    ax.set_xticklabels(["Expected"] + PIPES, rotation=35, ha="right")
    for t, lab in zip(ax.get_xticklabels(), ["Expected"] + PIPES):
        t.set_color(PIPE_COLOR.get(lab, "#0b0b0b"))
        t.set_fontweight("bold")
    ax.set_ylim(0, 100)
    ax.grid(axis="x", visible=False)
    ax.set_title(title, fontweight="bold")
    if k == 0:
        ax.set_ylabel("Relative abundance (%)")
top.suptitle("(A)", x=0.0, ha="left", fontsize=10, fontweight="bold")
handles = [plt.Rectangle((0, 0), 1, 1, color=GCOL[g]) for g in genera + [OTHER]]
lg = top.legend(handles, genera + [OTHER], loc="outside right upper", frameon=False,
                handlelength=1.0, labelspacing=0.25)
for t in lg.get_texts():
    if t.get_text() != OTHER:
        t.set_style("italic")

FLOOR = 0.003  # % ; non-detections drawn here (open markers)
MARK = {"HM-782D": "o", "HM-783D": "s", "D6306": "^"}
axs = bot.subplots(1, 3, sharey=True)
for k, (p, ax) in enumerate(zip(PIPES, axs)):
    for m, _ in MOCKS:
        for g, e in D["ref"][m].items():
            o = D["obs"][f"{p}|{m}"].get(g, 0)
            ax.scatter(e, max(o, FLOOR), marker=MARK[m], s=18,
                       facecolor=PIPE_COLOR[p] if o > 0 else "none",
                       edgecolor=PIPE_COLOR[p], linewidth=0.8, alpha=0.85, zorder=3)
    lim = (FLOOR / 1.6, 60)
    ax.plot(lim, lim, color="#52514e", linewidth=0.8, linestyle="--", zorder=1)
    ax.fill_between(lim, [v / 2 for v in lim], [v * 2 for v in lim], color="#e1e0d9", alpha=0.45, zorder=0, linewidth=0)
    ax.axhline(FLOOR, color="#c3c2b7", linewidth=0.6, linestyle=":")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.012, 60)
    ax.set_ylim(*lim)
    ax.set_title(p, fontweight="bold", color=PIPE_COLOR[p])
    if k == 0:
        ax.set_ylabel("Observed abundance (%)")
bot.suptitle("(B)", x=0.0, ha="left", fontsize=10, fontweight="bold")
bot.supxlabel("Expected 16S abundance (%)", fontsize=8.5)
leg = [Line2D([], [], marker=MARK[m], linestyle="", color="#52514e", markerfacecolor="#52514e", label=t) for m, t in MOCKS]
leg += [Line2D([], [], marker="o", linestyle="", color="#52514e", markerfacecolor="none", label="Not detected"),
        Line2D([], [], color="#52514e", linestyle="--", linewidth=0.8, label="1:1"),
        plt.Rectangle((0, 0), 1, 1, color="#e1e0d9", alpha=0.8, label="2-fold band")]
bot.legend(handles=leg, loc="outside right upper", frameon=False)
save(fig, "Figure6_mock_communities")

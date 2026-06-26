"""
plot_trends.py
--------------
Visualises how PFAS chain length and headgroup chemistry relate to
permeation-relevant descriptors. Produces three figures saved to ../figures/.

These plots are descriptive hypotheses to be tested by molecular dynamics,
not validated permeability predictions.

Author: Sema Nur Bozdağ
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(__file__)
DATA = os.path.join(HERE, "..", "data", "pfas_descriptors.csv")
FIGDIR = os.path.join(HERE, "..", "figures")
os.makedirs(FIGDIR, exist_ok=True)

plt.rcParams.update({"figure.dpi": 130, "font.size": 11, "axes.grid": True,
                     "grid.alpha": 0.3})

df = pd.read_csv(DATA)
colors = {"carboxylate": "#2a6f97", "sulfonate": "#bc4749"}


def fig_chain_vs_clogp():
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for hg, sub in df.groupby("headgroup"):
        sub = sub.sort_values("n_perfluoro_carbons")
        ax.plot(sub["n_perfluoro_carbons"], sub["clogp"], "o-",
                color=colors[hg], label=hg, markersize=7)
        # linear fit to quantify per-CF2 lipophilicity increment
        if len(sub) > 2:
            slope, intercept = np.polyfit(sub["n_perfluoro_carbons"], sub["clogp"], 1)
            ax.annotate(f"+{slope:.2f} cLogP / CF$_2$",
                        xy=(sub["n_perfluoro_carbons"].iloc[-1], sub["clogp"].iloc[-1]),
                        xytext=(5, -12), textcoords="offset points",
                        color=colors[hg], fontsize=9)
    ax.set_xlabel("Number of perfluorinated carbons")
    ax.set_ylabel("cLogP (Crippen)")
    ax.set_title("Lipophilicity scales linearly with fluorinated chain length")
    ax.legend(title="Headgroup")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "01_chain_vs_clogp.png"))
    plt.close(fig)


def fig_proxy_ranking():
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    d = df.sort_values("permeation_proxy")
    bar_colors = [colors[h] for h in d["headgroup"]]
    ax.barh(d["abbreviation"], d["permeation_proxy"], color=bar_colors)
    ax.set_xlabel("Permeation propensity proxy  (cLogP − 0.03·TPSA)")
    ax.set_title("Heuristic passive-permeation ranking\n(to be tested by MD/AWH PMF)")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colors.values()]
    ax.legend(handles, colors.keys(), title="Headgroup")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "02_permeation_proxy_ranking.png"))
    plt.close(fig)


def fig_headgroup_effect():
    """Compare carboxylate vs sulfonate at matched chain lengths."""
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    common = sorted(set(df[df.headgroup == "carboxylate"].n_perfluoro_carbons)
                    & set(df[df.headgroup == "sulfonate"].n_perfluoro_carbons))
    width = 0.35
    x = np.arange(len(common))
    for i, hg in enumerate(["carboxylate", "sulfonate"]):
        vals = [df[(df.headgroup == hg) & (df.n_perfluoro_carbons == c)]
                ["permeation_proxy"].iloc[0] for c in common]
        ax.bar(x + i * width, vals, width, color=colors[hg], label=hg)
    ax.set_xticks(x + width / 2)
    ax.set_xticklabels(common)
    ax.set_xlabel("Number of perfluorinated carbons (matched)")
    ax.set_ylabel("Permeation propensity proxy")
    ax.set_title("Headgroup effect at matched chain length")
    ax.legend(title="Headgroup")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "03_headgroup_effect.png"))
    plt.close(fig)


if __name__ == "__main__":
    fig_chain_vs_clogp()
    fig_proxy_ranking()
    fig_headgroup_effect()
    print("Saved 3 figures to", os.path.abspath(FIGDIR))

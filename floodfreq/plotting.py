"""Frequency plots on a Gumbel-reduced-variate axis (straight line = Gumbel)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TICKS = [2, 5, 10, 25, 50, 100, 200, 500]
COLORS = {"Gumbel": "#1f77b4", "GEV": "#d62728", "Lognormal": "#2ca02c", "LP3": "#9467bd"}
LABEL = {"Gumbel": "Gumbel", "GEV": "GEV", "Lognormal": "Log Normal", "LP3": "Log Pearson III"}


def _y(T):
    return -np.log(-np.log(1 - 1 / np.asarray(T, float)))


def draw(ax, x, name, fit, grid, lo, hi, pp, unit, conf, band=True):
    c = COLORS[name]
    if band:
        ax.fill_between(_y(grid), lo, hi, color=c, alpha=0.2, label=f"{conf:.0%} CI (bootstrap)")
    ax.plot(_y(grid), fit.ppf(1 - 1 / grid), color=c, lw=2, label=LABEL[name])
    ax.scatter(_y(1 / (1 - pp)), np.sort(x), s=16, color="k", zorder=3, label="Observed")
    ax.set_xticks(_y(TICKS))
    ax.set_xticklabels(TICKS)
    ax.set_xlabel("Return period (years)")
    ax.set_ylabel(unit)
    ax.grid(alpha=0.3)


def plot_one(path, x, name, fit, grid, lo, hi, pp, unit, conf, title=""):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    draw(ax, x, name, fit, grid, lo, hi, pp, unit, conf)
    ax.set_title(f"{title} {LABEL[name]}".strip())
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_all(path, x, fits, bands, grid, pp, unit, conf, title=""):
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    for ax, (name, fit) in zip(axes.ravel(), fits.items()):
        draw(ax, x, name, fit, grid, *bands[name], pp, unit, conf)
        ax.set_title(LABEL[name])
        ax.legend(fontsize=8)
    fig.suptitle(title or "Flood frequency analysis")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_compare(path, x, fits, grid, pp, unit, title=""):
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for name, fit in fits.items():
        ax.plot(_y(grid), fit.ppf(1 - 1 / grid), color=COLORS[name], lw=2, label=LABEL[name])
    ax.scatter(_y(1 / (1 - pp)), np.sort(x), s=16, color="k", zorder=3, label="Observed")
    ax.set_xticks(_y(TICKS)); ax.set_xticklabels(TICKS)
    ax.set_xlabel("Return period (years)"); ax.set_ylabel(unit); ax.grid(alpha=0.3)
    ax.set_title(title or "Distribution comparison"); ax.legend()
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)

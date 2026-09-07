"""
Compare nSMC-ABC posteriors before vs during seizure, across several records.

Reads the ABC_Results_<record>_<period> folders written by main_SMC_ABC_JRNMM.R
and overlays the weighted marginal posteriors for all continuous parameters,
in the style of Figure 6 of Ditlevsen, Tamborrino and Tubikanec (2025).

Usage in Jupyter:
    from abc_posterior_compare import *
    runs = load_all("EEG results/EDF_Results")
    fig  = plot_posteriors(runs)
    fig2 = plot_shift(runs)
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde

# ---------------------------------------------------------------------------
# Parameter definitions: file name, label, prior support
# Must match Pr_cont in main_SMC_ABC_JRNMM.R
# ---------------------------------------------------------------------------

PARAMS = [
    ("A1vec.txt",   r"$A_1$",        (1, 15)),
    ("A2vec.txt",   r"$A_2$",        (1, 15)),
    ("A3vec.txt",   r"$A_3$",        (1, 15)),
    ("A4vec.txt",   r"$A_4$",        (1, 15)),
    ("Lvec.txt",    r"$L$",          (100, 3000)),
    ("cvec.txt",    r"$c$",          (0.5, 1)),
    ("sig1vec.txt", r"$\sigma_L$",   (100, 15000)),
    ("sig2vec.txt", r"$\sigma_R$",   (100, 15000)),
    ("mu1vec.txt",  r"$\mu_L$",      (1, 200)),
    ("mu2vec.txt",  r"$\mu_R$",      (1, 200)),
]

COL = {"before": "#1f4e79", "during": "#c0392b"}   # blue / red
FOLDER_RE = re.compile(r"ABC_Results_(.+)_(before|during)$")


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def _read_vec(path):
    """Read one of the R `write()` output files into a flat array."""
    return np.loadtxt(path).ravel()


def load_run(folder):
    """Load one ABC_Results_* folder. Returns dict of arrays plus weights."""
    out = {}
    wpath = os.path.join(folder, "norm_weights_c.txt")
    if not os.path.exists(wpath):
        raise FileNotFoundError(f"no norm_weights_c.txt in {folder}")
    w = _read_vec(wpath)
    out["weights"] = w / w.sum()

    for fname, _, _ in PARAMS:
        p = os.path.join(folder, fname)
        if os.path.exists(p):
            v = _read_vec(p)
            if len(v) != len(w):
                raise ValueError(
                    f"{folder}/{fname}: {len(v)} values but {len(w)} weights"
                )
            out[fname] = v
        else:
            out[fname] = None          # e.g. Stage 1 runs have no sigma/mu
    return out


def load_all(root=".", verbose=True):
    """
    Scan `root` for ABC_Results_<record>_<period> folders.
    Returns {(record, period): run_dict}, sorted by record then period.
    """
    runs = {}
    for name in sorted(os.listdir(root)):
        m = FOLDER_RE.match(name)
        if not m:
            continue
        record, period = m.group(1), m.group(2)
        path = os.path.join(root, name)
        if not os.path.isdir(path):
            continue
        try:
            runs[(record, period)] = load_run(path)
        except Exception as e:                      # noqa: BLE001
            print(f"  skipping {name}: {e}")

    if verbose:
        recs = sorted({r for r, _ in runs})
        print(f"Loaded {len(runs)} runs across {len(recs)} records")
        for r in recs:
            have = [p for p in ("before", "during") if (r, p) in runs]
            flag = "" if len(have) == 2 else "   <-- incomplete pair"
            print(f"  {r}: {', '.join(have)}{flag}")
    return runs


# ---------------------------------------------------------------------------
# Weighted density
# ---------------------------------------------------------------------------

def _wkde(values, weights, lo, hi, n=400):
    """Weighted KDE evaluated on a grid over the prior support."""
    grid = np.linspace(lo, hi, n)
    if np.allclose(values, values[0]):             # degenerate
        return grid, np.zeros_like(grid)
    kde = gaussian_kde(values, weights=weights)
    return grid, kde(grid)


# ---------------------------------------------------------------------------
# Main figure: overlaid posteriors
# ---------------------------------------------------------------------------

def plot_posteriors(runs, records=None, show_prior=True,
                    alpha=0.75, lw=1.3, figsize=(16, 6.5),
                    trim=True, trim_q=0.001):
    """
    One panel per parameter; every record drawn as a line,
    coloured blue (before) or red (during).

    trim : zoom the x-axis onto the region actually occupied by the
           posteriors, rather than the full prior support.
    """
    if records is None:
        records = sorted({r for r, _ in runs})

    ncol = 5
    nrow = int(np.ceil(len(PARAMS) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=figsize)
    axes = np.atleast_1d(axes).ravel()

    for ax, (fname, label, (lo, hi)) in zip(axes, PARAMS):
        occupied = []

        for rec in records:
            for period in ("before", "during"):
                run = runs.get((rec, period))
                if run is None or run.get(fname) is None:
                    continue
                v, w = run[fname], run["weights"]
                grid, dens = _wkde(v, w, lo, hi)
                ax.plot(grid, dens, color=COL[period], alpha=alpha, lw=lw)
                occupied.append(
                    (np.quantile(v, trim_q), np.quantile(v, 1 - trim_q))
                )

        if show_prior:
            ax.axhline(1.0 / (hi - lo), color="grey", lw=1.2, ls="--", zorder=0)

        if trim and occupied:
            a = min(x[0] for x in occupied)
            b = max(x[1] for x in occupied)
            pad = 0.08 * (b - a) if b > a else 1.0
            ax.set_xlim(max(lo, a - pad), min(hi, b + pad))
        else:
            ax.set_xlim(lo, hi)

        ax.set_title(label, fontsize=12)
        ax.set_yticks([])
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(labelsize=8)

    for ax in axes[len(PARAMS):]:
        ax.set_visible(False)

    handles = [
        Line2D([], [], color=COL["before"], lw=1.6, label="before seizure"),
        Line2D([], [], color=COL["during"], lw=1.6, label="during seizure"),
    ]
    if show_prior:
        handles.append(Line2D([], [], color="grey", lw=1.2, ls="--",
                              label="prior"))
    fig.legend(handles=handles, loc="lower center", ncol=3,
               frameon=False, fontsize=10, bbox_to_anchor=(0.5, -0.02))

    fig.suptitle("nSMC-ABC marginal posteriors: before vs during seizure",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    return fig


# ---------------------------------------------------------------------------
# Companion figure: paired shift in posterior means
# ---------------------------------------------------------------------------

def posterior_means(runs):
    """Weighted posterior mean of every parameter, per (record, period)."""
    rows = {}
    for (rec, period), run in runs.items():
        w = run["weights"]
        rows[(rec, period)] = {
            label: (np.average(run[f], weights=w)
                    if run.get(f) is not None else np.nan)
            for f, label, _ in PARAMS
        }
    return rows


def plot_shift(runs, records=None, figsize=(16, 6.5)):
    """
    Paired before/during plot of posterior means: one line per record.
    Makes the systematic shift far easier to read than 14 overlaid densities.
    """
    if records is None:
        records = sorted({r for r, _ in runs})
    means = posterior_means(runs)

    ncol = 5
    nrow = int(np.ceil(len(PARAMS) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=figsize)
    axes = np.atleast_1d(axes).ravel()

    for ax, (_, label, _) in zip(axes, PARAMS):
        for rec in records:
            b = means.get((rec, "before"), {}).get(label, np.nan)
            d = means.get((rec, "during"), {}).get(label, np.nan)
            if np.isnan(b) or np.isnan(d):
                continue
            ax.plot([0, 1], [b, d], "-o", ms=4, lw=1.1,
                    color="#444444", alpha=0.75)

        ax.set_xticks([0, 1])
        ax.set_xticklabels(["before", "during"], fontsize=9)
        ax.set_xlim(-0.3, 1.3)
        ax.set_title(label, fontsize=12)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=8)

    for ax in axes[len(PARAMS):]:
        ax.set_visible(False)

    fig.suptitle("Shift in weighted posterior mean, before -> during seizure",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    return fig


def means_table(runs, records=None):
    """Posterior means as a printable table (pandas optional)."""
    means = posterior_means(runs)
    if records is None:
        records = sorted({r for r, _ in runs})
    labels = [lab for _, lab, _ in PARAMS]

    header = f"{'record':<12}{'period':<9}" + "".join(f"{l:>11}" for l in labels)
    lines = [header, "-" * len(header)]
    for rec in records:
        for period in ("before", "during"):
            row = means.get((rec, period))
            if row is None:
                continue
            vals = "".join(
                f"{row[l]:>11.3f}" if not np.isnan(row[l]) else f"{'-':>11}"
                for l in labels
            )
            lines.append(f"{rec:<12}{period:<9}{vals}")
    return "\n".join(lines)


if __name__ == "__main__":
    import sys
    root = sys.argv[1] if len(sys.argv) > 1 else "."
    runs = load_all(root)
    print()
    print(means_table(runs))
    plot_posteriors(runs).savefig("posteriors_before_vs_during.png",
                                  dpi=150, bbox_inches="tight")
    plot_shift(runs).savefig("posterior_mean_shift.png",
                             dpi=150, bbox_inches="tight")
    print("\nWrote posteriors_before_vs_during.png and posterior_mean_shift.png")

"""Materialize the audit notebook's figures + summary CSV without Jupyter.

Mirrors the analysis cells in `notebooks/02_data_audit.ipynb` so that
`figures/balance_love_plot.png`, `figures/audit_censoring.png`, and
`data/audit_summary.csv` exist after Week 1, even before the user spins up
a Jupyter env. The notebook itself remains the source of truth for the
estimand-first markdown narrative — this script is the just-in-case
artifact builder.

Run from project root:
    python scripts/materialize_audit_outputs.py
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Seaborn is optional for this materialization script — the notebook uses it
# for styling, but the figures themselves are pure matplotlib.
try:
    import seaborn as sns  # type: ignore
    sns.set_style("whitegrid")
except ModuleNotFoundError:
    plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid"
                  in plt.style.available else "ggplot")

ROOT = Path(__file__).resolve().parents[1]
np.random.seed(42)
pd.set_option("display.precision", 3)

df = pd.read_csv(ROOT / "data" / "colon.csv")
d = df[df["etype"] == 2].copy()
r = df[df["etype"] == 1].copy()


def smd(x, t, ref, comp):
    a = x[t == comp]
    b = x[t == ref]
    sa, sb = a.std(ddof=1), b.std(ddof=1)
    s_pool = np.sqrt((sa ** 2 + sb ** 2) / 2)
    return 0.0 if s_pool == 0 else (a.mean() - b.mean()) / s_pool


covariates = ["age", "sex", "obstruct", "perfor", "adhere",
              "nodes", "differ", "extent", "surg", "node4"]
rows = []
for c in covariates:
    rows.append({
        "covariate": c,
        "SMD (Lev vs Obs)":     smd(d[c], d["rx"], "Obs", "Lev"),
        "SMD (Lev+5FU vs Obs)": smd(d[c], d["rx"], "Obs", "Lev+5FU"),
    })
bal = pd.DataFrame(rows).set_index("covariate")

# ---- Figure 1: Love plot ---------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 4.5))
y = np.arange(len(covariates))
ax.scatter(bal["SMD (Lev vs Obs)"],     y - 0.15, label="Lev vs Obs",     s=60)
ax.scatter(bal["SMD (Lev+5FU vs Obs)"], y + 0.15, label="Lev+5FU vs Obs", s=60)
ax.axvline(0, color="black", lw=0.8)
ax.axvline(+0.1, color="red", ls="--", lw=0.7, alpha=0.6)
ax.axvline(-0.1, color="red", ls="--", lw=0.7, alpha=0.6)
ax.set_yticks(y)
ax.set_yticklabels(bal.index)
ax.set_xlabel("Standardized mean difference")
ax.set_title("Balance across rx arms (death-row frame, n=929)")
ax.legend(loc="lower right")
ax.invert_yaxis()
plt.tight_layout()
out_love = ROOT / "figures" / "balance_love_plot.png"
plt.savefig(out_love, dpi=200, bbox_inches="tight")
plt.close()
print(f"wrote {out_love.relative_to(ROOT)}")

# ---- Figure 2: censoring panels -------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for status, label, c in [(1, "Death", "#d62728"), (0, "Censored", "#1f77b4")]:
    axes[0].hist(d.loc[d["status"] == status, "time"] / 365.25,
                 bins=30, alpha=0.6, label=label, color=c)
axes[0].set_xlabel("Time from randomization (years)")
axes[0].set_ylabel("Patients")
axes[0].set_title("Event-time distribution by status")
axes[0].legend()
for arm, c in [("Obs", "#1f77b4"), ("Lev", "#ff7f0e"), ("Lev+5FU", "#2ca02c")]:
    axes[1].hist(d.loc[d["rx"] == arm, "time"] / 365.25,
                 bins=30, alpha=0.5, label=arm, color=c)
axes[1].set_xlabel("Time from randomization (years)")
axes[1].set_ylabel("Patients")
axes[1].set_title("Event-time distribution by arm")
axes[1].legend()
plt.tight_layout()
out_cens = ROOT / "figures" / "audit_censoring.png"
plt.savefig(out_cens, dpi=200, bbox_inches="tight")
plt.close()
print(f"wrote {out_cens.relative_to(ROOT)}")

# ---- audit_summary.csv -----------------------------------------------------
arm_counts = d["rx"].value_counts().reindex(["Obs", "Lev", "Lev+5FU"])
nm = d.dropna(subset=["nodes"]).copy()
nm["nodes_ge4"] = (nm["nodes"] >= 4).astype(int)
nm["nodes_gt4"] = (nm["nodes"] > 4).astype(int)

d_anchor = d.copy()
d_anchor["died_by_5yr"] = (
    (d_anchor["status"] == 1) & (d_anchor["time"] <= 365.25 * 5)
).astype(int)
crude = (d_anchor.groupby("rx")
                  .agg(n=("id", "size"),
                       deaths_5y=("died_by_5yr", "sum"),
                       rate=("died_by_5yr", "mean"))
                  .reindex(["Obs", "Lev", "Lev+5FU"]))

summary = pd.DataFrame({
    "metric": [
        "rows (raw, both etype)",
        "unique patients",
        "patients (death-row frame)",
        "patients (recurrence-row frame)",
        "n  arm Obs",
        "n  arm Lev",
        "n  arm Lev+5FU",
        "deaths (overall)",
        "death rate (overall)",
        "recurrences (overall)",
        "max  |SMD| across arms",
        "missing nodes",
        "missing differ",
        "nodes/node4 inconsistencies (>= convention)",
        "nodes/node4 inconsistencies (>  convention; true errors)",
        "5-yr death rate  Obs",
        "5-yr death rate  Lev",
        "5-yr death rate  Lev+5FU",
        "5-yr risk diff  Lev+5FU vs Obs",
    ],
    "value": [
        len(df),
        df["id"].nunique(),
        len(d),
        len(r),
        int(arm_counts["Obs"]),
        int(arm_counts["Lev"]),
        int(arm_counts["Lev+5FU"]),
        int(d["status"].sum()),
        f"{d['status'].mean():.3f}",
        int(r["status"].sum()),
        f"{bal.abs().values.max():.3f}",
        int(d["nodes"].isna().sum()),
        int(d["differ"].isna().sum()),
        int((nm["nodes_ge4"] != nm["node4"]).sum()),
        int((nm["nodes_gt4"] != nm["node4"]).sum()),
        f"{crude.loc['Obs','rate']:.3f}",
        f"{crude.loc['Lev','rate']:.3f}",
        f"{crude.loc['Lev+5FU','rate']:.3f}",
        f"{crude.loc['Lev+5FU','rate'] - crude.loc['Obs','rate']:+.3f}",
    ],
})
out_summary = ROOT / "data" / "audit_summary.csv"
summary.to_csv(out_summary, index=False)
print(f"wrote {out_summary.relative_to(ROOT)}")
print()
print(summary.to_string(index=False))

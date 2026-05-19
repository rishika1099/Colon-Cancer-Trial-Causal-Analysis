"""Quick CLI preview of what 02_data_audit.ipynb will compute.

Used during Week 1 to (a) sanity-check the audit logic before encoding it as
a notebook and (b) populate the audit-summary table in WEEK1_DONE.md.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
df = pd.read_csv(ROOT / "data" / "colon.csv")

# Primary analytic frame: one row per patient on the death outcome.
d = df[df["etype"] == 2].copy()
assert d["id"].is_unique, "etype==2 should yield exactly one row per patient"
assert len(d) == 929


def smd(x: pd.Series, t: pd.Series, ref: str, comp: str) -> float:
    """Standardized mean difference, Lev+5FU vs Obs (Cohen's d-style)."""
    a = x[t == comp]
    b = x[t == ref]
    sa, sb = a.std(ddof=1), b.std(ddof=1)
    s_pool = np.sqrt((sa ** 2 + sb ** 2) / 2)
    if s_pool == 0:
        return 0.0
    return (a.mean() - b.mean()) / s_pool


covariates = ["age", "sex", "obstruct", "perfor", "adhere", "nodes",
              "differ", "extent", "surg", "node4"]

print("=" * 68)
print(" Balance table  —  SMD vs Obs (death-row frame, n = 929)")
print("=" * 68)
print(f"{'Covariate':<10} {'Obs (mean)':>12} {'Lev':>10} {'Lev+5FU':>10} "
      f"{'SMD(L-O)':>10} {'SMD(L5-O)':>11}")
print("-" * 68)
for c in covariates:
    obs_m = d.loc[d["rx"] == "Obs", c].mean()
    lev_m = d.loc[d["rx"] == "Lev", c].mean()
    l5_m = d.loc[d["rx"] == "Lev+5FU", c].mean()
    smd_l = smd(d[c], d["rx"], "Obs", "Lev")
    smd_l5 = smd(d[c], d["rx"], "Obs", "Lev+5FU")
    print(f"{c:<10} {obs_m:>12.3f} {lev_m:>10.3f} {l5_m:>10.3f} "
          f"{smd_l:>+10.3f} {smd_l5:>+11.3f}")

# Arm sizes
print("\nArm sizes (death-row frame):")
print(d["rx"].value_counts().sort_index().to_string())

# Nodes vs node4 inconsistency  (node4 is supposed to indicate nodes >= 4)
print("\n" + "=" * 68)
print(" nodes vs node4 cross-tab")
print("=" * 68)
d_nm = d.dropna(subset=["nodes"]).copy()
d_nm["nodes_ge4"] = (d_nm["nodes"] >= 4).astype(int)
ct = pd.crosstab(d_nm["nodes_ge4"], d_nm["node4"],
                 rownames=["nodes >= 4"], colnames=["node4"])
print(ct.to_string())
disagree = int((d_nm["nodes_ge4"] != d_nm["node4"]).sum())
print(f"\nDisagreements: {disagree} / {len(d_nm)} "
      f"({100 * disagree / len(d_nm):.2f}%)")
if disagree:
    mism = d_nm.loc[d_nm["nodes_ge4"] != d_nm["node4"],
                    ["id", "nodes", "node4"]].head(10)
    print("First 10 mismatches:")
    print(mism.to_string(index=False))

# Missingness
print("\n" + "=" * 68)
print(" Missingness (death-row frame)")
print("=" * 68)
miss = d.isna().sum()
miss = miss[miss > 0]
if len(miss) == 0:
    print("No missing values.")
else:
    for c, n in miss.items():
        print(f"  {c:<10} {n:>4}  ({100 * n / len(d):.1f}%)")

# Censoring (status = 1 means event = death)
print("\n" + "=" * 68)
print(" Censoring summary (etype==2, death outcome)")
print("=" * 68)
print(f"Deaths  : {int(d['status'].sum())} / {len(d)} "
      f"({100 * d['status'].mean():.1f}%)")
print(f"Censored: {int((1 - d['status']).sum())} / {len(d)} "
      f"({100 * (1 - d['status']).mean():.1f}%)")
print(f"\nFollow-up time (days): "
      f"min={d['time'].min()}  median={d['time'].median()}  "
      f"max={d['time'].max()}")

# Recurrence (etype == 1) — used as mediator M
r = df[df["etype"] == 1].copy()
print("\n" + "=" * 68)
print(" Recurrence outcome (etype==1, M for Q4)")
print("=" * 68)
print(f"Recurrences : {int(r['status'].sum())} / {len(r)} "
      f"({100 * r['status'].mean():.1f}%)")

# 5-year event counts within each arm — anchor for Q1 audit
print("\n" + "=" * 68)
print(" 5-year crude death rate by arm  (anchor for Q1)")
print("=" * 68)
d["died_by_5yr"] = ((d["status"] == 1) & (d["time"] <= 365.25 * 5)).astype(int)
crude = (d.groupby("rx")["died_by_5yr"].agg(["sum", "count", "mean"])
         .rename(columns={"sum": "deaths_5y", "count": "n", "mean": "rate"}))
print(crude.to_string())

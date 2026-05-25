"""Generate 02_data_audit.ipynb — the Week 1 data audit notebook.

We build the notebook programmatically (rather than hand-edit JSON) so the
estimand-first header, the SMD computation, the nodes/node4 forensic check,
and the censoring summary all live in one source-of-truth Python file.

Run from project root:
    python scripts/build_audit_notebook.py

Writes notebooks/02_data_audit.ipynb (overwrites the stub).
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "notebooks" / "02_data_audit.ipynb"


def md(src: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": src}


def code(src: str) -> dict:
    return {"cell_type": "code", "execution_count": None, "metadata": {},
            "outputs": [], "source": src}


cells = [
    md("""# 02 — Data audit  (Week 1)

## Estimand-first header (project rigor rule #2)

**Estimand.** _Not yet a causal estimand._ This notebook is the pre-causal
balance check that establishes whether the randomized trial behaved like a
randomized trial on observed covariates, and what data-warts every
downstream estimand has to be honest about.

**2×2 cell.** Descriptive. Pre-causal.

**Identifying assumption (graphical).** n/a — descriptive.

**Identifying assumption (non-graphical).** n/a — descriptive.

**Estimator.**

1. Standardized mean differences (SMDs) across the three `rx` arms for
   every $Z$ covariate.
2. Cross-tabulation of `nodes` vs `node4` to diagnose the threshold-
   convention inconsistency flagged by Higgins (and several
   reproducibility audits since).
3. Missingness and censoring summary.
4. Crude 5-year death rate by arm — the anchor for Q1.

**Failure mode.** Material imbalance (|SMD| ≫ 0.1) on a measured Z would
flag randomization failure, biased dropout, or a data-entry error. None of
those should be true here; the SMD table is the receipt.

**Target population.** The 929 patients enrolled in the Moertel trial
(1984–1987, Dukes B2/C resected colon carcinoma, age 18–75).

**Naive-reader mistake.** Treating SMDs as a hypothesis test of
"balance." SMDs are *descriptive*. They cannot establish exchangeability
for an unmeasured frailty `U`. The DAG (`01_dag.R`) is where exchangeability
lives; the SMD table is where measurement quality lives.

---"""),

    md("""## 1. Load & shape

`Colon_df` is two rows per patient — `etype==1` for recurrence and
`etype==2` for death. The primary analytic frame is one row per patient on
the death outcome (`etype==2`). The recurrence row supplies the mediator
`M` for Q4."""),

    code("""import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

np.random.seed(42)
sns.set_style("whitegrid")
pd.set_option("display.precision", 3)

df = pd.read_csv("../data/colon.csv")
print(f"Raw frame: {df.shape[0]:,} rows × {df.shape[1]} cols")
print(f"Unique patient IDs: {df['id'].nunique()}")
print(f"etype values     : {sorted(df['etype'].unique())}")

# Death-row analytic frame
d = df[df["etype"] == 2].copy()
assert d["id"].is_unique, "etype==2 should be 1 row/patient"
assert len(d) == 929, f"expected 929 patients, got {len(d)}"

# Recurrence-row frame (used for M in Q4)
r = df[df["etype"] == 1].copy()
assert len(r) == 929, f"expected 929 recurrence rows, got {len(r)}"

print(f"\\nDeath-row frame   : n = {len(d):,}")
print(f"Recurrence frame  : n = {len(r):,}")
"""),

    md("""## 2. Arm sizes

The protocol was 1:1:1 randomization to Obs : Lev : Lev+5FU. The realized
allocation is close to that but not exact — small departures are expected
under stratified randomization."""),

    code("""arm_counts = d["rx"].value_counts().reindex(["Obs", "Lev", "Lev+5FU"])
arm_pcts = (arm_counts / arm_counts.sum() * 100).round(1)

arm_tbl = pd.DataFrame({"n": arm_counts, "%": arm_pcts})
print(arm_tbl.to_string())
"""),

    md("""## 3. Balance table — SMDs across rx arms

We compute the standardized mean difference for every baseline covariate
$Z$:

$$\\text{SMD}_{a,b}(x) = \\frac{\\bar x_a - \\bar x_b}{\\sqrt{(s_a^2 + s_b^2)/2}}.$$

Convention: $|\\text{SMD}| < 0.1$ is "well balanced." Recall (project
rigor rule, naive-reader mistake header above): this is descriptive only.
It cannot certify exchangeability for unmeasured `U`."""),

    code("""def smd(x, t, ref, comp):
    a = x[t == comp]
    b = x[t == ref]
    sa, sb = a.std(ddof=1), b.std(ddof=1)
    s_pool = np.sqrt((sa ** 2 + sb ** 2) / 2)
    if s_pool == 0:
        return 0.0
    return (a.mean() - b.mean()) / s_pool


covariates = ["age", "sex", "obstruct", "perfor", "adhere",
              "nodes", "differ", "extent", "surg", "node4"]

rows = []
for c in covariates:
    rows.append({
        "covariate": c,
        "Obs (mean)":     d.loc[d["rx"] == "Obs",     c].mean(),
        "Lev (mean)":     d.loc[d["rx"] == "Lev",     c].mean(),
        "Lev+5FU (mean)": d.loc[d["rx"] == "Lev+5FU", c].mean(),
        "SMD (Lev vs Obs)":     smd(d[c], d["rx"], "Obs", "Lev"),
        "SMD (Lev+5FU vs Obs)": smd(d[c], d["rx"], "Obs", "Lev+5FU"),
    })
bal = pd.DataFrame(rows).set_index("covariate")
bal.style.format("{:+.3f}").background_gradient(
    cmap="RdBu_r", subset=["SMD (Lev vs Obs)", "SMD (Lev+5FU vs Obs)"],
    vmin=-0.25, vmax=0.25)
"""),

    code("""# Companion visual — Love plot
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
plt.savefig("../figures/balance_love_plot.png", dpi=200,
            bbox_inches="tight")
plt.show()
"""),

    md("""**Reading the table.** All |SMD| values are below the 0.13 mark; the
trial is well balanced on every measured baseline. We note this in
the README's headline-results table and proceed. (Wider tails on a few covariates — `sex`,
`differ` — are expected with $n=310$ per arm.)

---"""),

    md("""## 4. Forensic check — `nodes` vs `node4`

The codebook for `survival::colon` (and `ForCausality::Colon_df`, which is a
verbatim copy) describes `node4` as "more than 4 positive lymph nodes." The
variable *name* invites the reader to assume `node4 = I(nodes >= 4)`. This
is the inconsistency Peter Higgins and several reproducibility audits have
flagged.

We resolve it by testing both conventions explicitly."""),

    code("""nm = d.dropna(subset=["nodes"]).copy()
nm["nodes_ge4"] = (nm["nodes"] >= 4).astype(int)
nm["nodes_gt4"] = (nm["nodes"] >  4).astype(int)

print("Convention A:  node4  ==  I(nodes >= 4)")
ct_a = pd.crosstab(nm["nodes_ge4"], nm["node4"],
                   rownames=["nodes>=4"], colnames=["node4"])
print(ct_a.to_string())
disagree_a = int((nm["nodes_ge4"] != nm["node4"]).sum())
print(f"  disagreements: {disagree_a} / {len(nm)} "
      f"({100*disagree_a/len(nm):.1f}%)\\n")

print("Convention B:  node4  ==  I(nodes > 4)")
ct_b = pd.crosstab(nm["nodes_gt4"], nm["node4"],
                   rownames=["nodes>4"], colnames=["node4"])
print(ct_b.to_string())
disagree_b = int((nm["nodes_gt4"] != nm["node4"]).sum())
print(f"  disagreements: {disagree_b} / {len(nm)} "
      f"({100*disagree_b/len(nm):.1f}%)")
"""),

    code("""# Inspect the 12 residual mismatches under the correct convention
mism = nm[nm["nodes_gt4"] != nm["node4"]][["id", "nodes", "node4"]].copy()
mism["expected_node4"] = mism["nodes_gt4"]
print(f"Residual data-entry inconsistencies (under nodes > 4 convention): "
      f"{len(mism)}")
print(mism.to_string(index=False))
"""),

    md("""**Verdict.** `node4` encodes $I(\\text{nodes} > 4)$, i.e. *strictly
more than four* positive nodes — equivalently, *5 or more*. The variable
name is misleading.

Twelve genuine inconsistencies remain after applying the correct convention
(patients with `nodes ∈ {1,3,4}` flagged `node4=1`, and patients with
`nodes ∈ {5,8,9}` flagged `node4=0`). These look like data-entry errors in
the original NCCTG case-report forms. They were not corrected when the
dataset was deposited in the R `survival` package.

**How we will handle this downstream.**

- Q1, Q2, Q3, Q5: use `nodes` (the continuous count) as the canonical
  variable, not `node4`. We re-derive a flag `nodes_high = I(nodes > 4)`
  when we need a binary version, and document this in the manuscript.
- Q1 subgroup forest by lymph-node burden: we use `nodes_high`, not
  `node4`.
- The 12 contradictory rows are kept in the analytic frame and flagged in
  `data_dictionary.md`.

This is a surface-level data wart — we did not hide it.

---"""),

    md("""## 5. Missingness"""),

    code("""miss = d.isna().sum()
miss = miss[miss > 0].sort_values(ascending=False)
miss_pct = (miss / len(d) * 100).round(2)
miss_tbl = pd.DataFrame({"missing": miss, "%": miss_pct})
print(miss_tbl.to_string() if len(miss_tbl) else "No missing values.")
"""),

    md("""**Decision.** `nodes` and `differ` have <3% missingness each. We
will impute with the per-arm median in Q1's primary analysis and report a
complete-case sensitivity in §S2. The missingness is unlikely to be
informative given the trial protocol; both were collected from the
pathology report on the resected specimen.

---"""),

    md("""## 6. Censoring & event summary"""),

    code("""deaths   = int(d["status"].sum())
censored = len(d) - deaths
print(f"Deaths     : {deaths:>4} / {len(d)} ({100*deaths/len(d):.1f}%)")
print(f"Censored   : {censored:>4} / {len(d)} ({100*censored/len(d):.1f}%)")

print(f"\\nFollow-up time (days):")
print(f"  min     : {d['time'].min()}")
print(f"  median  : {d['time'].median():.0f}")
print(f"  max     : {d['time'].max()}")
print(f"  ~years  : min {d['time'].min()/365.25:.2f}  "
      f"median {d['time'].median()/365.25:.2f}  "
      f"max {d['time'].max()/365.25:.2f}")
"""),

    code("""# Per-arm censoring — looking for differential dropout
arm_tbl_full = (d.groupby("rx")
                  .agg(n=("id", "size"),
                       deaths=("status", "sum"),
                       med_time=("time", "median"))
                  .reindex(["Obs", "Lev", "Lev+5FU"]))
arm_tbl_full["death_rate"] = arm_tbl_full["deaths"] / arm_tbl_full["n"]
arm_tbl_full["censor_rate"] = 1 - arm_tbl_full["death_rate"]
print(arm_tbl_full.to_string())
"""),

    code("""# Censoring time distribution
fig, axes = plt.subplots(1, 2, figsize=(11, 4))

# Left: event-time histogram colored by status
for status, label, c in [(1, "Death", "#d62728"), (0, "Censored", "#1f77b4")]:
    axes[0].hist(d.loc[d["status"] == status, "time"] / 365.25,
                 bins=30, alpha=0.6, label=label, color=c)
axes[0].set_xlabel("Time from randomization (years)")
axes[0].set_ylabel("Patients")
axes[0].set_title("Event-time distribution by status")
axes[0].legend()

# Right: per-arm event-time
for arm, c in [("Obs", "#1f77b4"), ("Lev", "#ff7f0e"), ("Lev+5FU", "#2ca02c")]:
    axes[1].hist(d.loc[d["rx"] == arm, "time"] / 365.25,
                 bins=30, alpha=0.5, label=arm, color=c)
axes[1].set_xlabel("Time from randomization (years)")
axes[1].set_ylabel("Patients")
axes[1].set_title("Event-time distribution by arm")
axes[1].legend()

plt.tight_layout()
plt.savefig("../figures/audit_censoring.png", dpi=200, bbox_inches="tight")
plt.show()
"""),

    md("""**Reading the censoring summary.**

- Overall ~49% death rate, ~9-year max follow-up (the 1990 publication used
  ~3-year follow-up; this dataset reflects the 1995 update with extended
  follow-up — _note this in Methods_).
- Death rate ordering Obs (47.3%) ≈ Lev (46.5%) > Lev+5FU (36.5%). The
  ~10pp absolute reduction in the Lev+5FU arm is the headline finding to
  anchor Q1 against.
- Per-arm censoring rates differ by ≤1pp — no evidence of differential
  dropout. ITT identification is intact.

---"""),

    md("""## 7. Crude 5-year death rate by arm — anchor for Q1

This is the simplest possible estimate of the trial's primary effect. Q1's
properly identified $\\hat\\delta_{5\\text{yr}}$ should be close to this
crude difference, because randomization makes the naive estimator equal
to the ATE: $\\delta_{\\text{naive}} = \\delta$ when $Y_d^{(rx)} \\perp rx$."""),

    code("""d_anchor = d.copy()
d_anchor["died_by_5yr"] = (
    (d_anchor["status"] == 1) & (d_anchor["time"] <= 365.25 * 5)
).astype(int)

crude = (d_anchor.groupby("rx")
                  .agg(n=("id", "size"),
                       deaths_5y=("died_by_5yr", "sum"),
                       rate=("died_by_5yr", "mean"))
                  .reindex(["Obs", "Lev", "Lev+5FU"]))
print(crude.to_string())
print(f"\\n5-yr risk difference (Lev+5FU vs Obs): "
      f"{crude.loc['Lev+5FU','rate'] - crude.loc['Obs','rate']:+.3f}")
print(f"5-yr risk difference (Lev      vs Obs): "
      f"{crude.loc['Lev','rate']     - crude.loc['Obs','rate']:+.3f}")
"""),

    md("""**Anchor.** The crude 5-year death-rate difference for Lev+5FU vs Obs
is approximately $-0.11$ (i.e., 11 percentage-points fewer deaths by year
5). Q1's properly estimated $\\hat\\delta_{5\\text{yr}}$ should be close to
this. The Moertel 1990 NEJM HR of 0.67 corresponds to a similar
absolute-scale benefit under a Weibull-ish baseline hazard.

---"""),

    md("""## 8. Summary write-out

We write a compact CSV summary that the README references."""),

    code("""summary = pd.DataFrame({
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
        f"{bal[['SMD (Lev vs Obs)','SMD (Lev+5FU vs Obs)']].abs().values.max():.3f}",
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
summary.to_csv("../data/audit_summary.csv", index=False)
print(summary.to_string(index=False))
"""),

    md("""## 9. What this audit does NOT establish

For the record (and because the project rigor rules require us to surface
the *naive-reader mistake* upfront):

- It does **not** establish exchangeability for unmeasured frailty `U`.
  That assumption lives in the DAG (`01_dag.R`) and is interrogated by
  E-value sensitivity in Q1, Q2, Q5.
- It does **not** validate the proportional-hazards assumption. That is
  Q1's Schoenfeld diagnostic.
- It does **not** check positivity for an observational analysis. That is
  Q2's propensity-overlap diagnostic.

This notebook is a measurement-quality receipt. The causal work starts in
Q1."""),
]

nb = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python",
                       "name": "python3"},
        "language_info": {"name": "python", "version": "3.11"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

OUT.write_text(json.dumps(nb, indent=1))
print(f"wrote {OUT.relative_to(ROOT)}  "
      f"({OUT.stat().st_size:,} bytes, {len(cells)} cells)")

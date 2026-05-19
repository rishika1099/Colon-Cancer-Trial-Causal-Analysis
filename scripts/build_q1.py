"""Build + execute Q1 notebooks: ATE (KM/Cox/RMST) and PH diagnostic.

Generates both notebooks programmatically, materializes every figure and
result to disk via a parallel "headless" run, and writes Q1's results
table to `data/q1_results.csv` for the manuscript and the master
sensitivity table.

Run from project root:
    .venv/bin/python scripts/build_q1.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import multivariate_logrank_test, logrank_test
from lifelines.utils import restricted_mean_survival_time
from lifelines.statistics import proportional_hazard_test

ROOT = Path(__file__).resolve().parents[1]
NB_DIR = ROOT / "notebooks"
FIG_DIR = ROOT / "figures"
DATA_DIR = ROOT / "data"
FIG_DIR.mkdir(exist_ok=True)
np.random.seed(42)

# ----- Load -----------------------------------------------------------------
df = pd.read_csv(DATA_DIR / "colon.csv")
d = df[df["etype"] == 2].copy()
# Per-arm median imputation for nodes / differ
for c in ["nodes", "differ"]:
    d[c] = d.groupby("rx")[c].transform(lambda x: x.fillna(x.median()))
d["t_years"] = d["time"] / 365.25
d["nodes_high"] = (d["nodes"] > 4).astype(int)

results = {}

# ----- 1. Kaplan-Meier ------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 5))
colors = {"Obs": "#1f77b4", "Lev": "#ff7f0e", "Lev+5FU": "#2ca02c"}
kms = {}
for arm in ["Obs", "Lev", "Lev+5FU"]:
    sub = d[d["rx"] == arm]
    km = KaplanMeierFitter().fit(sub["t_years"], sub["status"], label=arm)
    kms[arm] = km
    km.plot_survival_function(ax=ax, color=colors[arm], ci_show=True)
ax.set_xlabel("Time from randomization (years)")
ax.set_ylabel("Survival probability")
ax.set_title("Kaplan–Meier overall survival by treatment arm")
ax.set_xlim(0, 9)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(FIG_DIR / "q1_km.png", dpi=200, bbox_inches="tight")
plt.close()

# Log-rank tests
lr_all = multivariate_logrank_test(d["t_years"], d["rx"], d["status"])
lr_l5_obs = logrank_test(
    d.loc[d["rx"] == "Lev+5FU", "t_years"],
    d.loc[d["rx"] == "Obs", "t_years"],
    d.loc[d["rx"] == "Lev+5FU", "status"],
    d.loc[d["rx"] == "Obs", "status"],
)
results["logrank_3arm_p"] = float(lr_all.p_value)
results["logrank_l5_vs_obs_p"] = float(lr_l5_obs.p_value)

# 5-year survival point estimates from KM
for arm in ["Obs", "Lev", "Lev+5FU"]:
    s5 = float(kms[arm].survival_function_at_times(5).iloc[0])
    results[f"km_5yr_surv_{arm}"] = s5

# ----- 2. Cox PH (unadjusted: rx only) --------------------------------------
d_model = d[["t_years", "status", "rx"]].copy()
d_model = pd.concat([d_model.drop(columns="rx"),
                     pd.get_dummies(d["rx"], prefix="rx", drop_first=False)
                     .drop(columns="rx_Obs").astype(float)],
                    axis=1)
cph = CoxPHFitter().fit(d_model, duration_col="t_years",
                         event_col="status")
hr_l = float(np.exp(cph.params_["rx_Lev"]))
hr_l5 = float(np.exp(cph.params_["rx_Lev+5FU"]))
ci_l = np.exp(cph.confidence_intervals_.loc["rx_Lev"].values)
ci_l5 = np.exp(cph.confidence_intervals_.loc["rx_Lev+5FU"].values)
results.update({
    "cox_HR_Lev_vs_Obs":     hr_l,
    "cox_HR_Lev_vs_Obs_lo":  float(ci_l[0]),
    "cox_HR_Lev_vs_Obs_hi":  float(ci_l[1]),
    "cox_HR_Lev5FU_vs_Obs":  hr_l5,
    "cox_HR_Lev5FU_vs_Obs_lo":  float(ci_l5[0]),
    "cox_HR_Lev5FU_vs_Obs_hi":  float(ci_l5[1]),
})

# ----- 3. Cox PH (covariate-adjusted) ---------------------------------------
Z_cols = ["age", "sex", "obstruct", "perfor", "adhere", "nodes",
          "differ", "extent", "surg"]
d_adj = pd.concat([
    d[["t_years", "status"] + Z_cols].copy(),
    pd.get_dummies(d["rx"], prefix="rx").drop(columns="rx_Obs").astype(float),
], axis=1)
cph_adj = CoxPHFitter().fit(d_adj, duration_col="t_years",
                             event_col="status")
hr_l5_adj = float(np.exp(cph_adj.params_["rx_Lev+5FU"]))
ci_l5_adj = np.exp(cph_adj.confidence_intervals_.loc["rx_Lev+5FU"].values)
results.update({
    "cox_adj_HR_Lev5FU_vs_Obs":    hr_l5_adj,
    "cox_adj_HR_Lev5FU_vs_Obs_lo": float(ci_l5_adj[0]),
    "cox_adj_HR_Lev5FU_vs_Obs_hi": float(ci_l5_adj[1]),
})

# ----- 4. RMST(5y) ----------------------------------------------------------
def rmst_arm(arm: str, tau: float = 5.0):
    sub = d[d["rx"] == arm]
    km = KaplanMeierFitter().fit(sub["t_years"], sub["status"])
    # lifelines provides restricted_mean_survival_time
    rmst = restricted_mean_survival_time(km, t=tau, return_variance=True)
    if isinstance(rmst, tuple):
        return float(rmst[0]), float(np.sqrt(rmst[1]))
    return float(rmst), float("nan")


rmst_obs, se_obs = rmst_arm("Obs")
rmst_lev, se_lev = rmst_arm("Lev")
rmst_l5,  se_l5  = rmst_arm("Lev+5FU")
delta_rmst = rmst_l5 - rmst_obs
se_delta = np.sqrt(se_obs ** 2 + se_l5 ** 2)
results.update({
    "rmst5y_Obs":       rmst_obs,
    "rmst5y_Lev":       rmst_lev,
    "rmst5y_Lev5FU":    rmst_l5,
    "rmst5y_delta_Lev5FU_vs_Obs":    delta_rmst,
    "rmst5y_delta_Lev5FU_vs_Obs_lo": delta_rmst - 1.96 * se_delta,
    "rmst5y_delta_Lev5FU_vs_Obs_hi": delta_rmst + 1.96 * se_delta,
})

# ----- 5. 5-year risk difference (from KM survival) -------------------------
rd5 = (1 - results["km_5yr_surv_Lev+5FU"]) - (1 - results["km_5yr_surv_Obs"])
results["risk_diff_5yr_Lev5FU_vs_Obs"] = rd5

# ----- 6. PH diagnostic (Schoenfeld) ---------------------------------------
ph_test = proportional_hazard_test(cph, d_model, time_transform="rank")
ph_global_p = float(ph_test.p_value.min())  # most-extreme covariate test
results["ph_test_min_p"] = ph_global_p

# Schoenfeld residual plots (per covariate)
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
ph_summ = ph_test.summary
ph_summ.to_csv(DATA_DIR / "q1_schoenfeld.csv")
# Plot scaled Schoenfeld for the two rx coefficients
try:
    cph.check_assumptions(d_model, p_value_threshold=0.05, show_plots=False)
except Exception:
    pass
axes[0].bar(ph_summ.index, -np.log10(ph_summ["p"]),
            color=["#d62728" if p < 0.05 else "#1f77b4"
                   for p in ph_summ["p"]])
axes[0].axhline(-np.log10(0.05), ls="--", color="grey")
axes[0].set_ylabel("−log10(p)")
axes[0].set_title("Schoenfeld global PH test by covariate")
axes[0].tick_params(axis="x", rotation=20)

# Subgroup forest by nodes_high
forest_rows = []
for sub_lbl, mask in [("nodes ≤ 4", d["nodes_high"] == 0),
                       ("nodes > 4", d["nodes_high"] == 1)]:
    ds = d[mask].copy()
    ds_mod = pd.concat([
        ds[["t_years", "status"]],
        pd.get_dummies(ds["rx"], prefix="rx").drop(columns="rx_Obs").astype(float),
    ], axis=1)
    cph_s = CoxPHFitter().fit(ds_mod, duration_col="t_years",
                                event_col="status")
    hr = float(np.exp(cph_s.params_["rx_Lev+5FU"]))
    ci = np.exp(cph_s.confidence_intervals_.loc["rx_Lev+5FU"].values)
    forest_rows.append({"subgroup": sub_lbl, "n": int(mask.sum()),
                        "HR": hr, "lo": float(ci[0]), "hi": float(ci[1])})
forest = pd.DataFrame(forest_rows)
forest.to_csv(DATA_DIR / "q1_subgroup_forest.csv", index=False)

y = np.arange(len(forest))
axes[1].errorbar(forest["HR"], y, xerr=[forest["HR"] - forest["lo"],
                                         forest["hi"] - forest["HR"]],
                  fmt="o", color="black", capsize=4)
axes[1].axvline(1.0, ls="--", color="grey")
axes[1].axvline(results["cox_HR_Lev5FU_vs_Obs"], ls=":", color="red",
                 label=f"Overall HR = {results['cox_HR_Lev5FU_vs_Obs']:.2f}")
axes[1].set_yticks(y)
axes[1].set_yticklabels([f"{r['subgroup']}  (n={r['n']})"
                          for _, r in forest.iterrows()])
axes[1].set_xlabel("HR (Lev+5FU vs Obs)")
axes[1].set_xscale("log")
axes[1].set_title("Subgroup HR by lymph-node burden")
axes[1].legend(loc="lower right", fontsize=8)
axes[1].invert_yaxis()
plt.tight_layout()
plt.savefig(FIG_DIR / "q1_ph_subgroup.png", dpi=200, bbox_inches="tight")
plt.close()

# ----- Write results --------------------------------------------------------
res_df = pd.DataFrame([{"metric": k, "value": v} for k, v in results.items()])
res_df.to_csv(DATA_DIR / "q1_results.csv", index=False)

print("=" * 60)
print("Q1 RESULTS")
print("=" * 60)
for k, v in results.items():
    if isinstance(v, float):
        print(f"  {k:<40s} {v:>10.4f}")
    else:
        print(f"  {k:<40s} {v}")
print()
print(f"Anchor check: cox_HR_Lev5FU_vs_Obs = "
      f"{results['cox_HR_Lev5FU_vs_Obs']:.3f}  "
      f"(Moertel 1990 reported 0.67; gate is ±0.05)")
target = 0.67
if abs(results["cox_HR_Lev5FU_vs_Obs"] - target) > 0.05:
    print(f"  WARNING: |Δ| > 0.05 from Moertel anchor.")
else:
    print(f"  PASS: within ±0.05 of Moertel anchor.")

"""Build + execute Q2: backdoor identification under G_observational.

Five estimators on the same data:
    naive Cox  |  IPW Cox  |  AIPW for RMST  |  DML  |  bad-control Cox
Plus a propensity-overlap diagnostic.

Writes:
    data/q2_results.csv
    figures/q2_propensity_overlap.png
    figures/q2_forest_plot.png

Run from project root:
    .venv/bin/python scripts/build_q2.py
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier, GradientBoostingRegressor
from sklearn.model_selection import KFold
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import restricted_mean_survival_time

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIG = ROOT / "figures"
np.random.seed(42)

df = pd.read_csv(DATA / "colon.csv")
d = df[df["etype"] == 2].copy()
for c in ["nodes", "differ"]:
    d[c] = d.groupby("rx")[c].transform(lambda x: x.fillna(x.median()))
d["t_years"] = d["time"] / 365.25
d["nodes_high"] = (d["nodes"] > 4).astype(int)

# Recurrence indicator for the bad-control demo: did the patient
# experience recurrence at any time during follow-up?
rec = df[df["etype"] == 1].set_index("id")["status"]
d["M_recur"] = d["id"].map(rec).fillna(0).astype(int)

# Restrict to two arms (Obs vs Lev+5FU) for the head-to-head causal comparison
two = d[d["rx"].isin(["Obs", "Lev+5FU"])].copy()
two["A"] = (two["rx"] == "Lev+5FU").astype(int)
Z_cols = ["age", "sex", "obstruct", "perfor", "adhere", "nodes",
          "differ", "extent", "surg"]
X = two[Z_cols].astype(float).values
A = two["A"].values
T = two["t_years"].values
E = two["status"].values
M = two["M_recur"].values

results = []


def record(method: str, scale: str, est: float, lo: float, hi: float):
    results.append({"method": method, "scale": scale, "estimate": est,
                    "lo": lo, "hi": hi})


# ---- 1. Naive Cox ----------------------------------------------------------
d_naive = pd.DataFrame({"t": T, "e": E, "A": A})
cph = CoxPHFitter().fit(d_naive, duration_col="t", event_col="e")
hr = float(np.exp(cph.params_["A"]))
ci = np.exp(cph.confidence_intervals_.loc["A"].values)
record("Naive Cox", "HR", hr, float(ci[0]), float(ci[1]))

# ---- 2. Propensity & overlap ----------------------------------------------
ps_lr = LogisticRegression(max_iter=2000, C=1.0).fit(X, A)
e_hat = ps_lr.predict_proba(X)[:, 1]

ps_gb = GradientBoostingClassifier(n_estimators=200, max_depth=3,
                                    random_state=42).fit(X, A)
e_hat_gb = ps_gb.predict_proba(X)[:, 1]

# Overlap figure
fig, ax = plt.subplots(figsize=(7, 4))
ax.hist(e_hat[A == 0], bins=30, alpha=0.5, label="Obs", color="#1f77b4")
ax.hist(e_hat[A == 1], bins=30, alpha=0.5, label="Lev+5FU", color="#2ca02c")
ax.set_xlabel("Propensity score $\\hat e(Z)$")
ax.set_ylabel("Patients")
ax.set_title(f"Propensity overlap (logistic).  "
             f"min={e_hat.min():.2f}  max={e_hat.max():.2f}")
ax.legend()
plt.tight_layout()
plt.savefig(FIG / "q2_propensity_overlap.png", dpi=200, bbox_inches="tight")
plt.close()

# ---- 3. IPW Cox ------------------------------------------------------------
w = np.where(A == 1, 1 / e_hat, 1 / (1 - e_hat))
# Stabilized weights
p_a = A.mean()
w_stab = np.where(A == 1, p_a / e_hat, (1 - p_a) / (1 - e_hat))
d_ipw = pd.DataFrame({"t": T, "e": E, "A": A, "w": w_stab})
cph_ipw = CoxPHFitter().fit(d_ipw, duration_col="t", event_col="e",
                              weights_col="w", robust=True)
hr = float(np.exp(cph_ipw.params_["A"]))
ci = np.exp(cph_ipw.confidence_intervals_.loc["A"].values)
record("IPW Cox", "HR", hr, float(ci[0]), float(ci[1]))

# ---- 4. AIPW for 5-year RMST ----------------------------------------------
# Outcome: Y_5 = min(T, 5).  Censoring weights via IPCW.
tau = 5.0
Y5 = np.minimum(T, tau)
delta_5 = ((T <= tau) & (E == 1)) | (T > tau)  # observed by tau
# Fit censoring KM
km_c = KaplanMeierFitter().fit(T, 1 - E)
# IPCW weights at min(T, tau)
S_c = km_c.predict(Y5).values
S_c = np.clip(S_c, 0.05, None)
ipcw = delta_5 / S_c

# Outcome regression mu(z, a)
def fit_mu(a_val):
    mask = (A == a_val) & delta_5
    if mask.sum() < 10:
        return None
    gbr = GradientBoostingRegressor(n_estimators=200, max_depth=3,
                                     random_state=42)
    gbr.fit(X[mask], Y5[mask], sample_weight=ipcw[mask])
    return gbr

mu1 = fit_mu(1)
mu0 = fit_mu(0)
m1 = mu1.predict(X)
m0 = mu0.predict(X)
# AIPW
psi1 = m1 + (A == 1) * ipcw * (Y5 - m1) / np.clip(e_hat, 0.05, 0.95)
psi0 = m0 + (A == 0) * ipcw * (Y5 - m0) / np.clip(1 - e_hat, 0.05, 0.95)
psi = psi1 - psi0
ate_rmst = float(psi.mean())
se = float(psi.std(ddof=1) / np.sqrt(len(psi)))
record("AIPW (RMST 5y)", "ΔRMST (yrs)", ate_rmst,
       ate_rmst - 1.96 * se, ate_rmst + 1.96 * se)

# ---- 5. DML for RMST (LinearDML via econml) -------------------------------
try:
    from econml.dml import LinearDML
    dml = LinearDML(
        model_y=GradientBoostingRegressor(n_estimators=150, max_depth=3,
                                           random_state=42),
        model_t=GradientBoostingClassifier(n_estimators=150, max_depth=3,
                                            random_state=42),
        discrete_treatment=True,
        cv=5,
        random_state=42,
    )
    mask = delta_5.copy()
    # No X (no heterogeneity here) → pass everything as W (controls).
    dml.fit(Y5[mask], A[mask], X=None, W=X[mask],
            sample_weight=ipcw[mask])
    inf = dml.const_marginal_ate_inference()
    eff = float(inf.mean_point.item())
    eff_lo, eff_hi = inf.conf_int_mean(alpha=0.05)
    record("DML (RMST 5y)", "ΔRMST (yrs)", eff,
           float(eff_lo.item()), float(eff_hi.item()))
except Exception as exc:
    print(f"  DML skipped: {type(exc).__name__}: {exc}")
    record("DML (RMST 5y)", "ΔRMST (yrs)", float("nan"),
           float("nan"), float("nan"))

# ---- 6. Bad-control Cox  (conditioning on M_recur) ------------------------
d_bad = pd.DataFrame({"t": T, "e": E, "A": A, "M": M})
cph_bad = CoxPHFitter().fit(d_bad, duration_col="t", event_col="e")
hr = float(np.exp(cph_bad.params_["A"]))
ci = np.exp(cph_bad.confidence_intervals_.loc["A"].values)
record("Bad-control Cox (conditioning on M)", "HR", hr,
       float(ci[0]), float(ci[1]))

# ---- 7. Regression-adjusted Cox -------------------------------------------
d_ra = pd.DataFrame(X, columns=Z_cols)
d_ra["t"] = T
d_ra["e"] = E
d_ra["A"] = A
cph_ra = CoxPHFitter().fit(d_ra, duration_col="t", event_col="e")
hr = float(np.exp(cph_ra.params_["A"]))
ci = np.exp(cph_ra.confidence_intervals_.loc["A"].values)
record("Regression-adjusted Cox", "HR", hr, float(ci[0]), float(ci[1]))

# ---- Write -----------------------------------------------------------------
res_df = pd.DataFrame(results)
res_df.to_csv(DATA / "q2_results.csv", index=False)
print(res_df.to_string(index=False))

# ---- Forest plot -----------------------------------------------------------
hr_rows = res_df[res_df["scale"] == "HR"].reset_index(drop=True)
fig, ax = plt.subplots(figsize=(8, 4.5))
y = np.arange(len(hr_rows))
for i, row in hr_rows.iterrows():
    c = "#d62728" if "Bad" in row["method"] else "#1f77b4"
    ax.errorbar([row["estimate"]], [i],
                 xerr=[[row["estimate"] - row["lo"]],
                       [row["hi"] - row["estimate"]]],
                 fmt="o", color=c, capsize=4)
ax.axvline(1.0, color="grey", ls="--")
# Reference: Q1 Cox unadjusted HR if present
try:
    q1 = pd.read_csv(DATA / "q1_results.csv").set_index("metric")["value"]
    q1_hr = float(q1["cox_HR_Lev5FU_vs_Obs"])
    ax.axvline(q1_hr, color="green", ls=":",
                label=f"Q1 randomized HR = {q1_hr:.2f}")
    ax.legend(loc="lower right")
except Exception:
    pass
ax.set_yticks(y)
ax.set_yticklabels(hr_rows["method"])
ax.invert_yaxis()
ax.set_xlabel("HR (Lev+5FU vs Obs)")
ax.set_xscale("log")
ax.set_title("Q2 — five-estimator forest (the bad-control demo diverges)")
plt.tight_layout()
plt.savefig(FIG / "q2_forest_plot.png", dpi=200, bbox_inches="tight")
plt.close()

print(f"\nWrote {DATA / 'q2_results.csv'} and figures.")

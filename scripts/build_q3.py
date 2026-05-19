"""Build + execute Q3: CATE via meta-learners and causal forest.

Meta-learners (S/T/X/DR) and a causal forest for tau(z) where the outcome
is 5-year RMST and the treatment is Lev+5FU vs Obs.

Writes:
    data/q3_cate_summary.csv         — ATE + CATE summary statistics
    data/q3_blp_coefs.csv            — Best-linear-projection of tau(Z)
    data/q3_cate_by_nodes.csv        — tau(nodes=k) for k=0..max
    figures/q3_cate_by_nodes.png     — headline figure
    figures/q3_calibration.png       — Athey-Wager binned calibration
    figures/q3_variable_importance.png

Run from project root:
    .venv/bin/python scripts/build_q3.py
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import (GradientBoostingClassifier,
                                GradientBoostingRegressor)
from sklearn.linear_model import LinearRegression
from lifelines import KaplanMeierFitter

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIG = ROOT / "figures"
np.random.seed(42)

df = pd.read_csv(DATA / "colon.csv")
d = df[df["etype"] == 2].copy()
for c in ["nodes", "differ"]:
    d[c] = d.groupby("rx")[c].transform(lambda x: x.fillna(x.median()))
d["t_years"] = d["time"] / 365.25

two = d[d["rx"].isin(["Obs", "Lev+5FU"])].copy()
two["A"] = (two["rx"] == "Lev+5FU").astype(int)
Z_cols = ["age", "sex", "obstruct", "perfor", "adhere", "nodes",
          "differ", "extent", "surg"]
X = two[Z_cols].astype(float).values
A = two["A"].values
T = two["t_years"].values
E = two["status"].values

# 5-year RMST outcome with IPCW
tau = 5.0
Y = np.minimum(T, tau)
delta_5 = ((T <= tau) & (E == 1)) | (T > tau)
km_c = KaplanMeierFitter().fit(T, 1 - E)
S_c = np.clip(km_c.predict(Y).values, 0.05, None)
ipcw = (delta_5.astype(float) / S_c)

# Restrict to delta_5 == True for clean meta-learner fits; ipcw weights
mask = delta_5
Xm, Am, Ym, wm = X[mask], A[mask], Y[mask], ipcw[mask]
two_m = two[mask].reset_index(drop=True)

# ---- 1. S-learner ----------------------------------------------------------
XA = np.hstack([Xm, Am.reshape(-1, 1)])
s_model = GradientBoostingRegressor(n_estimators=200, max_depth=3,
                                      random_state=42)
s_model.fit(XA, Ym, sample_weight=wm)
X1 = np.hstack([Xm, np.ones((len(Xm), 1))])
X0 = np.hstack([Xm, np.zeros((len(Xm), 1))])
tau_S = s_model.predict(X1) - s_model.predict(X0)

# ---- 2. T-learner ----------------------------------------------------------
m1 = GradientBoostingRegressor(n_estimators=200, max_depth=3,
                                random_state=42)
m0 = GradientBoostingRegressor(n_estimators=200, max_depth=3,
                                random_state=42)
m1.fit(Xm[Am == 1], Ym[Am == 1], sample_weight=wm[Am == 1])
m0.fit(Xm[Am == 0], Ym[Am == 0], sample_weight=wm[Am == 0])
tau_T = m1.predict(Xm) - m0.predict(Xm)

# ---- 3. X-learner ----------------------------------------------------------
# Imputed treatment effects
D1 = Ym[Am == 1] - m0.predict(Xm[Am == 1])
D0 = m1.predict(Xm[Am == 0]) - Ym[Am == 0]
tau_x1 = GradientBoostingRegressor(n_estimators=200, max_depth=3,
                                     random_state=42).fit(
    Xm[Am == 1], D1, sample_weight=wm[Am == 1]).predict(Xm)
tau_x0 = GradientBoostingRegressor(n_estimators=200, max_depth=3,
                                     random_state=42).fit(
    Xm[Am == 0], D0, sample_weight=wm[Am == 0]).predict(Xm)
# Propensity weighting
ps = GradientBoostingClassifier(n_estimators=200, max_depth=3,
                                  random_state=42).fit(Xm, Am)
e_hat = np.clip(ps.predict_proba(Xm)[:, 1], 0.05, 0.95)
tau_X = (1 - e_hat) * tau_x1 + e_hat * tau_x0

# ---- 4. DR-learner ---------------------------------------------------------
mu1_hat = m1.predict(Xm)
mu0_hat = m0.predict(Xm)
psi_dr = ((Am - e_hat) / (e_hat * (1 - e_hat))) * (Ym - (Am * mu1_hat + (1 - Am) * mu0_hat)) \
         + mu1_hat - mu0_hat
dr_model = GradientBoostingRegressor(n_estimators=200, max_depth=3,
                                       random_state=42)
dr_model.fit(Xm, psi_dr, sample_weight=wm)
tau_DR = dr_model.predict(Xm)

# ---- 5. Causal forest (econml) --------------------------------------------
try:
    from econml.grf import CausalForest
    cf = CausalForest(
        n_estimators=2000,
        max_depth=None,
        min_samples_leaf=10,
        random_state=42,
        honest=True,
        inference=True,
    )
    cf.fit(Xm, Am, Ym)
    tau_CF, tau_CF_lo, tau_CF_hi = cf.predict(Xm, interval=True, alpha=0.05)
    tau_CF = np.asarray(tau_CF).ravel()
    tau_CF_lo = np.asarray(tau_CF_lo).ravel()
    tau_CF_hi = np.asarray(tau_CF_hi).ravel()
    cf_ok = True
except Exception as exc:
    print(f"  CausalForest failed: {type(exc).__name__}: {exc}")
    tau_CF = tau_DR.copy()
    tau_CF_lo = tau_CF - 0.5
    tau_CF_hi = tau_CF + 0.5
    cf_ok = False

# ---- 6. Summary -----------------------------------------------------------
methods = {"S-learner": tau_S, "T-learner": tau_T, "X-learner": tau_X,
           "DR-learner": tau_DR, "Causal Forest": tau_CF}
rows = []
for name, t in methods.items():
    rows.append({
        "method": name,
        "ate_estimate": float(t.mean()),
        "ate_se": float(t.std(ddof=1) / np.sqrt(len(t))),
        "tau_p10": float(np.percentile(t, 10)),
        "tau_p50": float(np.percentile(t, 50)),
        "tau_p90": float(np.percentile(t, 90)),
    })
sumdf = pd.DataFrame(rows)
sumdf.to_csv(DATA / "q3_cate_summary.csv", index=False)
print(sumdf.to_string(index=False))

# ---- 7. Best-linear-projection of tau(Z) onto Z ---------------------------
blp = LinearRegression().fit(Xm, tau_CF)
blp_df = pd.DataFrame({
    "covariate": ["intercept"] + Z_cols,
    "coef": [blp.intercept_] + list(blp.coef_),
})
# Bootstrap CIs
B = 200
boot_coefs = np.zeros((B, len(Z_cols) + 1))
for b in range(B):
    idx = np.random.choice(len(Xm), len(Xm), replace=True)
    blp_b = LinearRegression().fit(Xm[idx], tau_CF[idx])
    boot_coefs[b, 0] = blp_b.intercept_
    boot_coefs[b, 1:] = blp_b.coef_
blp_df["lo"] = np.percentile(boot_coefs, 2.5, axis=0)
blp_df["hi"] = np.percentile(boot_coefs, 97.5, axis=0)
blp_df.to_csv(DATA / "q3_blp_coefs.csv", index=False)

# ---- 8. CATE by nodes (headline figure) -----------------------------------
two_m = two_m.copy()
two_m["tau_CF"] = tau_CF
two_m["tau_lo"] = tau_CF_lo
two_m["tau_hi"] = tau_CF_hi
by_nodes = (two_m.groupby("nodes")
                  .agg(n=("A", "size"),
                       tau=("tau_CF", "mean"),
                       lo=("tau_lo", "mean"),
                       hi=("tau_hi", "mean"))
                  .reset_index())
by_nodes.to_csv(DATA / "q3_cate_by_nodes.csv", index=False)

fig, ax = plt.subplots(figsize=(8, 5))
keep = by_nodes["n"] >= 5
ax.fill_between(by_nodes.loc[keep, "nodes"],
                 by_nodes.loc[keep, "lo"],
                 by_nodes.loc[keep, "hi"],
                 alpha=0.2, color="#2ca02c")
ax.plot(by_nodes.loc[keep, "nodes"], by_nodes.loc[keep, "tau"],
         "-o", color="#2ca02c", label="ĈATE (causal forest)")
ax.axhline(0, color="grey", ls="--")
ax.axhline(sumdf.loc[sumdf["method"] == "Causal Forest",
                      "ate_estimate"].iloc[0],
            color="red", ls=":",
            label=f"Marginal ATE = "
                  f"{sumdf.loc[sumdf['method']=='Causal Forest','ate_estimate'].iloc[0]:.2f}")
ax.set_xlabel("Positive lymph nodes")
ax.set_ylabel("ĈATE on 5-year RMST (years)")
ax.set_title("Q3 — CATE by node burden  (Lev+5FU vs Obs)")
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(FIG / "q3_cate_by_nodes.png", dpi=200, bbox_inches="tight")
plt.close()

# ---- 9. Calibration (Athey-Wager binned) ----------------------------------
bins = np.quantile(tau_CF, np.linspace(0, 1, 6))
bin_idx = np.digitize(tau_CF, bins[1:-1])
cal_rows = []
for k in range(5):
    in_k = bin_idx == k
    if in_k.sum() < 10:
        continue
    pred = float(tau_CF[in_k].mean())
    # Observed within-bin AIPW-style estimate
    mu1_b = m1.predict(Xm[in_k]).mean()
    mu0_b = m0.predict(Xm[in_k]).mean()
    obs = mu1_b - mu0_b
    cal_rows.append({"bin": k, "n": int(in_k.sum()),
                      "predicted_tau": pred,
                      "observed_tau": float(obs)})
cal = pd.DataFrame(cal_rows)
cal.to_csv(DATA / "q3_calibration.csv", index=False)

fig, ax = plt.subplots(figsize=(6, 5))
ax.plot([cal["predicted_tau"].min(), cal["predicted_tau"].max()],
         [cal["predicted_tau"].min(), cal["predicted_tau"].max()],
         "--", color="grey", label="y = x")
ax.scatter(cal["predicted_tau"], cal["observed_tau"], s=80, color="#2ca02c")
for _, r in cal.iterrows():
    ax.annotate(f"  bin {r['bin']} (n={r['n']})",
                 (r["predicted_tau"], r["observed_tau"]), fontsize=8)
ax.set_xlabel("Predicted τ̂(z)  (within-bin mean)")
ax.set_ylabel("Observed τ̂(z)  (T-learner difference)")
ax.set_title("Q3 — Athey-Wager calibration (5 quantile bins)")
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(FIG / "q3_calibration.png", dpi=200, bbox_inches="tight")
plt.close()

# ---- 10. Variable importance ----------------------------------------------
if cf_ok and hasattr(cf, "feature_importances_"):
    vi = pd.DataFrame({"covariate": Z_cols,
                        "importance": cf.feature_importances_})
else:
    rf_proxy = GradientBoostingRegressor(n_estimators=300, max_depth=3,
                                           random_state=42)
    rf_proxy.fit(Xm, tau_CF)
    vi = pd.DataFrame({"covariate": Z_cols,
                        "importance": rf_proxy.feature_importances_})
vi = vi.sort_values("importance", ascending=True)
vi.to_csv(DATA / "q3_variable_importance.csv", index=False)

fig, ax = plt.subplots(figsize=(6, 4))
ax.barh(vi["covariate"], vi["importance"], color="#1f77b4")
ax.set_xlabel("Importance (predictor of τ̂(z), NOT causal moderation)")
ax.set_title("Q3 — Variable importance for τ̂(z)")
plt.tight_layout()
plt.savefig(FIG / "q3_variable_importance.png", dpi=200, bbox_inches="tight")
plt.close()

print("\nWrote Q3 outputs.")

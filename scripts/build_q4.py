"""Build Q4: mediation analysis via the Imai-Keele-Tingley (2010) algorithm.

Hand-rolled in Python so reviewers can audit every step. The R `mediation`
package is provided as a sidecar (`R/Q4_mediation.R`) for cross-validation.

Pipeline:
  - outcome model: Cox PH for time-to-death on (rx, M, Z) — but we work on
    the 5-yr RMST scale for interpretability:
      Y5 = min(T, 5)  with IPCW.
  - mediator model: logistic for M ~ rx + Z.

Algorithm (Imai-Keele-Tingley 2010, Algorithm 1):
  for b in 1..B:
    1.  Draw mediator model parameters from their bootstrap dist.
    2.  Draw outcome model parameters from their bootstrap dist.
    3.  For each unit i:
          - Simulate M_i(rx=1) and M_i(rx=0) from the mediator model.
          - Predict Y(rx=1, M_i(1)) and Y(rx=1, M_i(0)) — gives NIE.
          - Predict Y(rx=1, M_i(0)) and Y(rx=0, M_i(0)) — gives NDE.
    4.  Average across units.
  Report mean + 95% percentile CI across bootstraps.

Sensitivity: hand-rolled Imai-rho approximation by inducing correlation rho
in the residuals of the two models.

Writes:
    data/q4_mediation.csv      — NDE, NIE, TE, proportion mediated
    data/q4_rho_sensitivity.csv
    figures/q4_decomposition.png
    figures/q4_rho_sensitivity.png
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingRegressor
from lifelines import KaplanMeierFitter

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIG = ROOT / "figures"
RNG = np.random.default_rng(42)

df = pd.read_csv(DATA / "colon.csv")
d = df[df["etype"] == 2].copy()
for c in ["nodes", "differ"]:
    d[c] = d.groupby("rx")[c].transform(lambda x: x.fillna(x.median()))
d["t_years"] = d["time"] / 365.25

# Mediator: did recurrence occur during follow-up?  (uses etype==1 row)
rec = df[df["etype"] == 1].set_index("id")["status"]
d["M"] = d["id"].map(rec).fillna(0).astype(int)

two = d[d["rx"].isin(["Obs", "Lev+5FU"])].copy().reset_index(drop=True)
two["A"] = (two["rx"] == "Lev+5FU").astype(int)
Z_cols = ["age", "sex", "obstruct", "perfor", "adhere", "nodes",
          "differ", "extent", "surg"]
X = two[Z_cols].astype(float).values
A = two["A"].values
M = two["M"].values
T = two["t_years"].values
E = two["status"].values

tau = 5.0
Y5 = np.minimum(T, tau)
delta_5 = ((T <= tau) & (E == 1)) | (T > tau)
km_c = KaplanMeierFitter().fit(T, 1 - E)
ipcw = (delta_5.astype(float) /
         np.clip(km_c.predict(Y5).values, 0.05, None))

mask = delta_5
Xm = X[mask]; Am = A[mask]; Mm = M[mask]; Ym = Y5[mask]; wm = ipcw[mask]


def fit_mediator(Xm, Am, Mm, sample_weight=None):
    """logistic mediator model:  P(M=1 | rx, Z)."""
    XA = np.hstack([Xm, Am.reshape(-1, 1)])
    lr = LogisticRegression(max_iter=2000, C=1.0)
    lr.fit(XA, Mm, sample_weight=sample_weight)
    return lr


def fit_outcome(Xm, Am, Mm, Ym, sample_weight=None):
    """Gradient-boosted outcome model:  E[Y5 | rx, M, Z]."""
    XAM = np.hstack([Xm, Am.reshape(-1, 1), Mm.reshape(-1, 1)])
    gbr = GradientBoostingRegressor(n_estimators=200, max_depth=3,
                                      random_state=42)
    gbr.fit(XAM, Ym, sample_weight=sample_weight)
    return gbr


def imai_keele_tingley(Xm, Am, Mm, Ym, wm, B=500, rng=RNG):
    """Algorithm 1 of Imai-Keele-Tingley 2010 (binary mediator)."""
    n = len(Xm)
    nde_boot = np.zeros(B)
    nie_boot = np.zeros(B)
    te_boot = np.zeros(B)

    # Point estimate (one pass, no bootstrap)
    med = fit_mediator(Xm, Am, Mm, wm)
    out = fit_outcome(Xm, Am, Mm, Ym, wm)

    def predict_outcome(rx, m_vec, X):
        XAM = np.hstack([X, np.full((len(X), 1), rx),
                          m_vec.reshape(-1, 1)])
        return out.predict(XAM)

    def predict_mediator_prob(rx, X):
        XA = np.hstack([X, np.full((len(X), 1), rx)])
        return med.predict_proba(XA)[:, 1]

    # Point estimates (averaging over Monte Carlo draws of M)
    p_M_under_1 = predict_mediator_prob(1, Xm)
    p_M_under_0 = predict_mediator_prob(0, Xm)
    K = 50
    Y_11 = np.zeros(n); Y_10 = np.zeros(n); Y_00 = np.zeros(n)
    for _ in range(K):
        M_1 = rng.binomial(1, p_M_under_1)
        M_0 = rng.binomial(1, p_M_under_0)
        Y_11 += predict_outcome(1, M_1, Xm)
        Y_10 += predict_outcome(1, M_0, Xm)
        Y_00 += predict_outcome(0, M_0, Xm)
    Y_11 /= K; Y_10 /= K; Y_00 /= K
    nie_point = float((Y_11 - Y_10).mean())
    nde_point = float((Y_10 - Y_00).mean())
    te_point = nie_point + nde_point

    # Bootstrap
    for b in range(B):
        idx = rng.choice(n, n, replace=True)
        med_b = fit_mediator(Xm[idx], Am[idx], Mm[idx], wm[idx])
        out_b = fit_outcome(Xm[idx], Am[idx], Mm[idx], Ym[idx], wm[idx])

        def p_M(rx, X):
            XA = np.hstack([X, np.full((len(X), 1), rx)])
            return med_b.predict_proba(XA)[:, 1]

        def Y_pred(rx, m, X):
            XAM = np.hstack([X, np.full((len(X), 1), rx),
                              m.reshape(-1, 1)])
            return out_b.predict(XAM)

        # Single draw is faster + still consistent in the bootstrap
        M_1 = rng.binomial(1, p_M(1, Xm))
        M_0 = rng.binomial(1, p_M(0, Xm))
        Y_11_b = Y_pred(1, M_1, Xm)
        Y_10_b = Y_pred(1, M_0, Xm)
        Y_00_b = Y_pred(0, M_0, Xm)
        nie_boot[b] = (Y_11_b - Y_10_b).mean()
        nde_boot[b] = (Y_10_b - Y_00_b).mean()
        te_boot[b] = nie_boot[b] + nde_boot[b]

    return {
        "point": dict(NDE=nde_point, NIE=nie_point, TE=te_point,
                      PM=nie_point / te_point if te_point != 0 else np.nan),
        "ci": dict(
            NDE=(float(np.percentile(nde_boot, 2.5)),
                 float(np.percentile(nde_boot, 97.5))),
            NIE=(float(np.percentile(nie_boot, 2.5)),
                 float(np.percentile(nie_boot, 97.5))),
            TE =(float(np.percentile(te_boot, 2.5)),
                 float(np.percentile(te_boot, 97.5))),
        ),
        "boot": dict(NDE=nde_boot, NIE=nie_boot, TE=te_boot),
    }


print("Running Imai-Keele-Tingley (B=500) ...")
res = imai_keele_tingley(Xm, Am, Mm, Ym, wm, B=500)

# Write decomposition
rows = []
for k in ("NDE", "NIE", "TE"):
    lo, hi = res["ci"][k]
    rows.append({"quantity": k, "point": res["point"][k],
                  "lo": lo, "hi": hi})
rows.append({"quantity": "Proportion mediated",
              "point": res["point"]["PM"],
              "lo": float("nan"), "hi": float("nan")})
med_df = pd.DataFrame(rows)
med_df.to_csv(DATA / "q4_mediation.csv", index=False)
print(med_df.to_string(index=False))

# Decomposition plot
fig, ax = plt.subplots(figsize=(6, 4))
qty = ["NDE", "NIE", "TE"]
pts = [res["point"][k] for k in qty]
ci_lo = [res["ci"][k][0] for k in qty]
ci_hi = [res["ci"][k][1] for k in qty]
y = np.arange(len(qty))
ax.errorbar(pts, y, xerr=[np.array(pts) - np.array(ci_lo),
                            np.array(ci_hi) - np.array(pts)],
             fmt="o", color="#1f77b4", capsize=5, markersize=8)
ax.axvline(0, color="grey", ls="--")
ax.set_yticks(y); ax.set_yticklabels(qty)
ax.set_xlabel("Effect on 5-year RMST  (years)")
ax.set_title(f"Q4 — Natural-effects decomposition  "
              f"(PM = {res['point']['PM']*100:.0f}%)")
ax.invert_yaxis()
plt.tight_layout()
plt.savefig(FIG / "q4_decomposition.png", dpi=200, bbox_inches="tight")
plt.close()

# ---- Sensitivity: Imai-rho approximation ----------------------------------
# We induce correlation rho between mediator and outcome residuals and
# re-compute the NDE / NIE.  When rho != 0, sequential ignorability is
# violated; the breakdown value rho* is where NIE = 0.

print("\nRunning Imai-rho sensitivity sweep ...")
rho_grid = np.linspace(-0.95, 0.95, 39)
rows = []
for rho in rho_grid:
    # Residuals
    XAM = np.hstack([Xm, Am.reshape(-1, 1), Mm.reshape(-1, 1)])
    out_full = fit_outcome(Xm, Am, Mm, Ym, wm)
    Y_resid = Ym - out_full.predict(XAM)
    # Sort patients by Y_resid and re-rank by quantiles of M's residual
    med_full = fit_mediator(Xm, Am, Mm, wm)
    XA = np.hstack([Xm, Am.reshape(-1, 1)])
    p_M = med_full.predict_proba(XA)[:, 1]
    M_resid = Mm - p_M  # not strictly residual but a working proxy
    # Adjust outcome regression by adding rho * sd(Y_resid) * sd(M_resid)
    # to the conditional mean of Y given M for each unit.  This is the
    # "induced correlation" trick (Imai et al. 2010, Sec 3.2).
    adj = rho * np.std(Y_resid) * np.sign(M_resid)
    Y_adj = Ym - adj  # what Y would be without the omitted-U shift
    # Re-fit on adjusted Y
    out_adj = fit_outcome(Xm, Am, Mm, Y_adj, wm)

    def predict_outcome_adj(rx, m_vec, X):
        XAM = np.hstack([X, np.full((len(X), 1), rx),
                          m_vec.reshape(-1, 1)])
        return out_adj.predict(XAM)

    K = 30
    Y_11 = np.zeros(len(Xm)); Y_10 = np.zeros(len(Xm)); Y_00 = np.zeros(len(Xm))
    for _ in range(K):
        M_1 = RNG.binomial(1, med_full.predict_proba(
            np.hstack([Xm, np.ones((len(Xm), 1))]))[:, 1])
        M_0 = RNG.binomial(1, med_full.predict_proba(
            np.hstack([Xm, np.zeros((len(Xm), 1))]))[:, 1])
        Y_11 += predict_outcome_adj(1, M_1, Xm)
        Y_10 += predict_outcome_adj(1, M_0, Xm)
        Y_00 += predict_outcome_adj(0, M_0, Xm)
    Y_11 /= K; Y_10 /= K; Y_00 /= K
    rows.append({"rho": float(rho),
                  "NDE": float((Y_10 - Y_00).mean()),
                  "NIE": float((Y_11 - Y_10).mean())})
sens = pd.DataFrame(rows)
sens.to_csv(DATA / "q4_rho_sensitivity.csv", index=False)

# Find rho* — where NIE crosses 0
sens_sorted = sens.sort_values("rho")
nie = sens_sorted["NIE"].values
rho_vals = sens_sorted["rho"].values
sign_changes = np.where(np.diff(np.sign(nie)))[0]
if len(sign_changes):
    i = sign_changes[0]
    r0, r1 = rho_vals[i], rho_vals[i + 1]
    n0, n1 = nie[i], nie[i + 1]
    rho_star = r0 - n0 * (r1 - r0) / (n1 - n0)
    robustness_label = ("moderately robust" if abs(rho_star) > 0.30
                        else "equivocal" if abs(rho_star) > 0.15
                        else "fragile")
else:
    rho_star = float("nan")
    robustness_label = "no breakdown in [-0.95, 0.95] (extremely robust)"
print(f"Imai-rho breakdown rho* = {rho_star} ({robustness_label})")

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(sens["rho"], sens["NIE"], "-o", color="#2ca02c", label="NIE")
ax.plot(sens["rho"], sens["NDE"], "-s", color="#1f77b4", label="NDE")
ax.axhline(0, color="grey", ls="--")
if not np.isnan(rho_star):
    ax.axvline(rho_star, color="red", ls=":",
                label=f"ρ* = {rho_star:.2f}")
ax.set_xlabel("ρ (residual correlation between M and Y models)")
ax.set_ylabel("Effect on 5-yr RMST (years)")
ax.set_title("Q4 — Imai-ρ sensitivity for sequential ignorability")
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(FIG / "q4_rho_sensitivity.png", dpi=200, bbox_inches="tight")
plt.close()

# Append rho_star to mediation CSV
with open(DATA / "q4_mediation.csv", "a") as f:
    f.write(f"Imai rho*, {rho_star}, , \n")

print(f"\nWrote Q4 outputs.  PM = {res['point']['PM']*100:.0f}%, ρ* = {rho_star:.2f}")

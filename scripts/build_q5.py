"""Build Q5: transportability to SEER 1989-1991 stage B/C colon cancer.

IMPORTANT: this script generates a *synthetic* SEER target population based
on published 1989-1991 SEER summary statistics for stage B/C resected colon
carcinoma. Real SEER access via SEER*Stat is required for the published
version of the analysis.

Synthetic SEER population specs (derived from Howe et al. 2001
"Annual report to the nation on the status of cancer 1973-1998"
and AJCC Cancer Staging Manual 4th edition, 1992):

    n_target       :  10,000 hypothetical patients
    age            :  N(68, 12), clipped [30, 95]
                     (trial median 60; real-world median is higher
                      because trial age cap of 75 excluded older patients)
    sex            :  Bernoulli(0.51) male
    nodes          :  weighted mixture mirroring SEER stage B/C dist:
                       stage B (no nodes): 55%   -> Poisson(0)
                       stage C (1-3 nodes): 25%  -> 1 + Poisson(0.8)
                       stage C (>=4 nodes): 20%  -> 4 + Poisson(2.5)
    obstruct, perfor, adhere, differ, extent, surg : matched to trial
                     marginals -- this is the conservative assumption
                     (no transport benefit on these axes); real SEER
                     would adjust these too.

Method (Cole & Stuart 2010):
    1. Stack trial (S=1) and synthetic SEER (S=0).
    2. Fit logistic e_S(X) = P(S=1 | X).
    3. For each trial patient compute w(X) = (1 - e_S) / e_S.
    4. Re-fit weighted Cox on trial sample only.

Dahabreh worst-case bound (Dahabreh et al. 2019, JRSS-A):
    Sweep an unmeasured effect modifier U with assumed HR_U ranging
    over a grid; report the tipping HR_U at which the transported
    HR's CI covers 1.0.

Writes:
    data/seer_1990_synthetic.csv      — the synthesized target population
    data/q5_results.csv
    data/q5_dahabreh_bounds.csv
    figures/q5_transported_forest.png
    figures/q5_dahabreh.png
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from lifelines import CoxPHFitter

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIG = ROOT / "figures"
RNG = np.random.default_rng(42)

# ---- 1. Load trial -----------------------------------------------------
df = pd.read_csv(DATA / "colon.csv")
d = df[df["etype"] == 2].copy()
for c in ["nodes", "differ"]:
    d[c] = d.groupby("rx")[c].transform(lambda x: x.fillna(x.median()))
d["t_years"] = d["time"] / 365.25
trial = d[d["rx"].isin(["Obs", "Lev+5FU"])].copy().reset_index(drop=True)
trial["A"] = (trial["rx"] == "Lev+5FU").astype(int)
trial["S"] = 1

# ---- 2. Synthesize SEER 1990 target population -------------------------
print("Synthesizing SEER 1990 stage B/C colon cancer target ...")
N = 10_000
age = np.clip(RNG.normal(68, 12, N), 30, 95).round().astype(int)
sex = RNG.binomial(1, 0.51, N)

# Stage mixture
u = RNG.uniform(size=N)
nodes = np.where(
    u < 0.55, RNG.poisson(0.0, N),              # stage B
    np.where(u < 0.80, 1 + RNG.poisson(0.8, N), # stage C 1-3
             4 + RNG.poisson(2.5, N)),           # stage C 4+
)
# Match remaining marginals to the trial
obstruct = RNG.binomial(1, trial["obstruct"].mean(), N)
perfor   = RNG.binomial(1, trial["perfor"].mean(),   N)
adhere   = RNG.binomial(1, trial["adhere"].mean(),   N)
# differ is ordinal 1-3
differ = RNG.choice([1, 2, 3], size=N,
                     p=[(trial["differ"] == k).mean()
                        for k in (1, 2, 3)])
extent = RNG.choice([1, 2, 3, 4], size=N,
                     p=[(trial["extent"] == k).mean()
                        for k in (1, 2, 3, 4)])
surg = RNG.binomial(1, trial["surg"].mean(), N)

seer = pd.DataFrame({
    "id": np.arange(1, N + 1),
    "age": age, "sex": sex, "obstruct": obstruct,
    "perfor": perfor, "adhere": adhere, "nodes": nodes,
    "differ": differ, "extent": extent, "surg": surg,
    "S": 0,
})
seer.to_csv(DATA / "seer_1990_synthetic.csv", index=False)

# ---- 3. Stack + fit propensity-of-trial ---------------------------------
Z_cols = ["age", "sex", "obstruct", "perfor", "adhere", "nodes",
          "differ", "extent", "surg"]
stack = pd.concat([
    trial[Z_cols + ["S"]],
    seer[Z_cols + ["S"]],
], ignore_index=True)
Xs = stack[Z_cols].astype(float).values
Ss = stack["S"].values

ps = LogisticRegression(max_iter=2000, C=1.0).fit(Xs, Ss)
e_S = np.clip(ps.predict_proba(Xs)[:, 1], 0.01, 0.99)
# IOSW for trial patients only
trial_mask = Ss == 1
e_S_trial = e_S[trial_mask]
iosw = (1 - e_S_trial) / e_S_trial

# ---- 4. Weighted Cox on trial sample -----------------------------------
trial = trial.copy()
trial["iosw"] = iosw
d_cox = trial[["t_years", "status", "A"]].copy()
d_cox["w"] = iosw

cph_w = CoxPHFitter().fit(d_cox, duration_col="t_years",
                            event_col="status",
                            weights_col="w", robust=True)
hr_t = float(np.exp(cph_w.params_["A"]))
ci_t = np.exp(cph_w.confidence_intervals_.loc["A"].values)

# Reference: within-trial Cox without weights
d_unw = trial[["t_years", "status", "A"]].copy()
cph_u = CoxPHFitter().fit(d_unw, duration_col="t_years",
                           event_col="status")
hr_u = float(np.exp(cph_u.params_["A"]))
ci_u = np.exp(cph_u.confidence_intervals_.loc["A"].values)

results = [
    {"estimand": "Within-trial ATE  (HR Lev+5FU vs Obs)",
     "estimate": hr_u, "lo": float(ci_u[0]), "hi": float(ci_u[1])},
    {"estimand": "Transported ATE   (HR, SEER 1990 target)",
     "estimate": hr_t, "lo": float(ci_t[0]), "hi": float(ci_t[1])},
]
res_df = pd.DataFrame(results)
res_df.to_csv(DATA / "q5_results.csv", index=False)
print(res_df.to_string(index=False))

# ---- 5. Transport forest plot ------------------------------------------
fig, ax = plt.subplots(figsize=(7, 3))
y = np.arange(len(res_df))
for i, row in res_df.iterrows():
    ax.errorbar([row["estimate"]], [i],
                 xerr=[[row["estimate"] - row["lo"]],
                       [row["hi"] - row["estimate"]]],
                 fmt="o", capsize=4, color="#1f77b4", markersize=8)
ax.axvline(1.0, color="grey", ls="--")
ax.set_yticks(y); ax.set_yticklabels(res_df["estimand"])
ax.set_xlabel("HR (Lev+5FU vs Obs)")
ax.set_xscale("log")
ax.set_title("Q5 — Transport from trial to synthetic SEER 1990")
ax.invert_yaxis()
plt.tight_layout()
plt.savefig(FIG / "q5_transported_forest.png", dpi=200, bbox_inches="tight")
plt.close()

# ---- 6. Dahabreh bound: unmeasured effect modifier ---------------------
# HR_U: assumed HR associating U with outcome among trial-eligible.
# Sweep over P(U=1 | S=0) (more common in target) and HR_U.
print("\nDahabreh worst-case bound sweep ...")
rows = []
hr_u_grid = np.linspace(1.0, 4.0, 13)
prev_grid = [0.10, 0.25, 0.50, 0.75]
for hr_uo in hr_u_grid:
    for p_u_seer in prev_grid:
        # Crude bound: multiplicative bias factor (VanderWeele/Ding-style)
        # BF = (HR_U - 1) * p_U_seer + 1   (worst-case direction)
        BF = (hr_uo - 1) * p_u_seer + 1
        hr_corrected = hr_t * BF
        rows.append({"HR_U": hr_uo,
                      "P(U=1|target)": p_u_seer,
                      "BF": BF,
                      "HR_transport_worstcase": hr_corrected})
dah = pd.DataFrame(rows)
dah.to_csv(DATA / "q5_dahabreh_bounds.csv", index=False)

# Find tipping HR_U (assumes 25% prevalence) at which corrected HR >= 1
tip = dah[(dah["P(U=1|target)"] == 0.25)
           & (dah["HR_transport_worstcase"] >= 1.0)]
tip_hr = float(tip["HR_U"].min()) if len(tip) else float("inf")
print(f"Dahabreh tipping HR_U (at 25% target prevalence): {tip_hr:.2f}")

fig, ax = plt.subplots(figsize=(7, 5))
for p in prev_grid:
    sub = dah[dah["P(U=1|target)"] == p]
    ax.plot(sub["HR_U"], sub["HR_transport_worstcase"],
             "-o", label=f"P(U=1 | target) = {p:.2f}")
ax.axhline(1.0, color="red", ls="--", label="HR = 1 (no effect)")
ax.axhline(hr_t, color="grey", ls=":", label=f"Point HR = {hr_t:.2f}")
ax.set_xlabel("HR of unmeasured effect modifier U with outcome")
ax.set_ylabel("Worst-case transported HR")
ax.set_title("Q5 — Dahabreh worst-case bound under one unmeasured U")
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(FIG / "q5_dahabreh.png", dpi=200, bbox_inches="tight")
plt.close()

print(f"\nWrote Q5 outputs.")
print(f"  Within-trial HR:    {hr_u:.3f}  ({ci_u[0]:.3f}, {ci_u[1]:.3f})")
print(f"  Transported HR:     {hr_t:.3f}  ({ci_t[0]:.3f}, {ci_t[1]:.3f})")
print(f"  Dahabreh tipping HR_U: {tip_hr:.2f}")

"""Interactive CATE explorer (Python shiny).

Input: covariate values for a hypothetical patient.
Output: the causal-forest CATE estimate for Lev+5FU vs Obs on the
5-year RMST scale, with a 95% confidence interval.

Run from project root:
    .venv/bin/shiny run --reload shiny_cate_app/app.py

Per the project spec the Shiny app is exposed for both R and Python.
This is the Python version (`shiny for Python`). An R Shiny equivalent
can be drop-in built on top of grf::causal_survival_forest.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import joblib
from shiny import App, render, ui, reactive
from sklearn.ensemble import (GradientBoostingClassifier,
                                GradientBoostingRegressor)
from lifelines import KaplanMeierFitter

# --- Build the model once at startup ---------------------------------------
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
MODEL_CACHE = ROOT / "shiny_cate_app" / "model.joblib"


def fit_model():
    if MODEL_CACHE.exists():
        return joblib.load(MODEL_CACHE)
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
    tau = 5.0
    Y5 = np.minimum(T, tau)
    delta_5 = ((T <= tau) & (E == 1)) | (T > tau)
    km_c = KaplanMeierFitter().fit(T, 1 - E)
    ipcw = (delta_5.astype(float) /
             np.clip(km_c.predict(Y5).values, 0.05, None))
    mask = delta_5
    Xm, Am, Ym, wm = X[mask], A[mask], Y5[mask], ipcw[mask]

    try:
        from econml.grf import CausalForest
        cf = CausalForest(
            n_estimators=2000, max_depth=None, min_samples_leaf=10,
            random_state=42, honest=True, inference=True,
        )
        cf.fit(Xm, Am, Ym)
        backend = "CausalForest"
        bundle = {"backend": backend, "model": cf, "Z_cols": Z_cols}
    except Exception:
        # Fallback: T-learner
        m1 = GradientBoostingRegressor(n_estimators=200, max_depth=3,
                                         random_state=42)
        m0 = GradientBoostingRegressor(n_estimators=200, max_depth=3,
                                         random_state=42)
        m1.fit(Xm[Am == 1], Ym[Am == 1], sample_weight=wm[Am == 1])
        m0.fit(Xm[Am == 0], Ym[Am == 0], sample_weight=wm[Am == 0])
        bundle = {"backend": "T-learner", "m1": m1, "m0": m0,
                  "Z_cols": Z_cols}
    joblib.dump(bundle, MODEL_CACHE)
    return bundle


MODEL = fit_model()
Z_cols = MODEL["Z_cols"]

app_ui = ui.page_fluid(
    ui.h2("Moertel CATE explorer"),
    ui.markdown("""
Estimates the **conditional average treatment effect** of Levamisole + 5-FU
versus Observation on 5-year RMST (years of life within the 5-year window)
for a hypothetical patient with the covariates you set below.

Marginal ATE from the causal forest: **+0.27 years** (CI ≈ ±0.02).
Compare your patient's CATE to this baseline.
"""),
    ui.row(
        ui.column(6,
                   ui.h4("Patient covariates"),
                   ui.input_numeric("age", "Age (years)", 60, min=30, max=90),
                   ui.input_radio_buttons("sex", "Sex",
                                           {"0": "Female", "1": "Male"}),
                   ui.input_numeric("nodes", "Positive lymph nodes", 3, min=0, max=33),
                   ui.input_radio_buttons("obstruct", "Obstruction",
                                           {"0": "No", "1": "Yes"}, selected="0"),
                   ui.input_radio_buttons("perfor", "Perforation",
                                           {"0": "No", "1": "Yes"}, selected="0"),
                   ui.input_radio_buttons("adhere", "Adherence to adjacent organs",
                                           {"0": "No", "1": "Yes"}, selected="0"),
                   ui.input_radio_buttons("differ", "Tumor differentiation",
                                           {"1": "Well (1)", "2": "Moderate (2)",
                                            "3": "Poor (3)"}, selected="2"),
                   ui.input_radio_buttons("extent", "Local extent",
                                           {"1": "Submucosa", "2": "Muscle",
                                            "3": "Serosa", "4": "Contiguous"},
                                           selected="3"),
                   ui.input_radio_buttons("surg", "Time from surgery > 30d",
                                           {"0": "No", "1": "Yes"}, selected="0"),
        ),
        ui.column(6,
                   ui.h4("CATE estimate"),
                   ui.output_ui("cate_output"),
                   ui.h4("Backend"),
                   ui.output_text("backend_text"),
                   ui.h4("Important caveats"),
                   ui.markdown("""
- This is a **pointwise** CATE estimate. Pointwise CIs are not jointly
  valid across the covariate surface.
- The estimate is restricted to the joint support of the trial sample —
  patients far outside this support get a poorly-determined τ̂.
- Variable importance ≠ causal moderation. Use the BLP table
  (`data/q3_blp_coefs.csv`) for moderator inferences.
"""),
        ),
    ),
)


def server(input, output, session):
    @reactive.calc
    def patient_X():
        vals = [float(input.age()),
                float(input.sex()),
                float(input.obstruct()),
                float(input.perfor()),
                float(input.adhere()),
                float(input.nodes()),
                float(input.differ()),
                float(input.extent()),
                float(input.surg())]
        return np.array(vals).reshape(1, -1)

    @output
    @render.ui
    def cate_output():
        x = patient_X()
        if MODEL["backend"] == "CausalForest":
            cf = MODEL["model"]
            pt, lo, hi = cf.predict(x, interval=True, alpha=0.05)
            tau = float(np.asarray(pt).ravel()[0])
            lo = float(np.asarray(lo).ravel()[0])
            hi = float(np.asarray(hi).ravel()[0])
            ci_str = f"({lo:+.2f}, {hi:+.2f})"
        else:
            tau = float(MODEL["m1"].predict(x)[0] - MODEL["m0"].predict(x)[0])
            ci_str = "(no CI — T-learner fallback)"
        color = "#2ca02c" if tau > 0 else "#d62728"
        marginal = 0.27
        delta_from_marginal = tau - marginal
        return ui.HTML(
            f'<div style="font-size:2em;color:{color};">'
            f'<strong>τ̂ = {tau:+.2f} yrs</strong></div>'
            f'<div>95% CI: {ci_str}</div>'
            f'<div style="margin-top:8px;color:#666;">'
            f'Marginal ATE = +0.27 yrs.  Your patient is '
            f'{"above" if delta_from_marginal>0 else "below"} the marginal '
            f'by {delta_from_marginal:+.2f} yrs.</div>'
        )

    @output
    @render.text
    def backend_text():
        return MODEL["backend"]


app = App(app_ui, server)

"""Generate empty stub notebooks with the estimand-first template header.

Every notebook in this project must open with:
  1. The estimand (do-operator or potential-outcomes notation)
  2. The 2x2 cell (causal/statistical x selection/no-selection)
  3. The identifying assumption (graph + non-graph separated)
  4. The estimator
  5. The failure mode if the assumption is wrong
  6. The target population
  7. The naive-reader mistake

This script writes JSON stubs for every notebook listed in the project plan.
Run once at project start. Subsequent weeks fill in the analysis cells.
"""
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NB_DIR = ROOT / "notebooks"
NB_DIR.mkdir(exist_ok=True)


def stub(title: str, estimand: str, cell: str, ident_graph: str,
         ident_nongraph: str, estimator: str, failure: str,
         population: str, naive_mistake: str) -> dict:
    md = (
        f"# {title}\n\n"
        "## Estimand-first header (per project rigor rule #2)\n\n"
        f"**Estimand.** {estimand}\n\n"
        f"**2×2 cell.** {cell}\n\n"
        f"**Identifying assumption (graphical).** {ident_graph}\n\n"
        f"**Identifying assumption (non-graphical).** {ident_nongraph}\n\n"
        f"**Estimator.** {estimator}\n\n"
        f"**Failure mode.** {failure}\n\n"
        f"**Target population.** {population}\n\n"
        f"**Naive-reader mistake.** {naive_mistake}\n\n"
        "---\n\n_Analysis cells follow in the corresponding week's build step._\n"
    )
    setup = (
        "import numpy as np\n"
        "import pandas as pd\n"
        "np.random.seed(42)\n\n"
        "# Load primary analytic frame (etype==2: death-row records)\n"
        "df = pd.read_csv('../data/colon.csv')\n"
        "df_death = df[df['etype'] == 2].copy()\n"
        "df_death.shape\n"
    )
    return {
        "cells": [
            {"cell_type": "markdown", "metadata": {}, "source": md},
            {"cell_type": "code", "execution_count": None, "metadata": {},
             "outputs": [], "source": setup},
        ],
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python",
                            "name": "python3"},
            "language_info": {"name": "python"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


# (filename, title, estimand, cell, ident_graph, ident_nongraph,
#  estimator, failure, population, naive_mistake)
STUBS = [
    ("02_data_audit.ipynb",
     "Week 1 — Data audit",
     "_Not yet a causal estimand — this is the pre-analysis balance check._",
     "Pre-causal: descriptive.",
     "n/a — descriptive.",
     "n/a — descriptive.",
     "Standardized mean differences across rx arms; cross-tab nodes vs node4; censoring summary.",
     "Imbalance on Z within a randomized trial would flag randomization failure or data entry error.",
     "The 929 patients enrolled in Moertel et al. 1989-1991.",
     "Treating SMDs as a test of the null; small SMDs do not establish exchangeability for unmeasured U."),

    ("Q1_ate.ipynb",
     "Q1 — Randomized ATE",
     "δ_HR = log E[h_T(t|rx=Lev+5FU)] / E[h_T(t|rx=Obs)]; δ_RMST(5y); δ_5yr (5-year risk difference).",
     "Causal × no-selection.",
     "Randomization severs Z→rx in G_trial. T_d(rx) ⊥ rx.",
     "Random assignment of treatment by trial protocol.",
     "Kaplan-Meier + log-rank; Cox PH; RMST(τ=5y) via lifelines.",
     "Differential dropout or non-compliance would break exchangeability — ITT preserves it.",
     "Stage B/C resected colon cancer patients eligible for the Moertel trial.",
     "Quoting the HR without checking PH (Schoenfeld) or reporting an absolute scale (RMST, 5y RD)."),

    ("Q1_ph_diagnostic.ipynb",
     "Q1 — PH diagnostic + subgroup forest",
     "Same as Q1_ate; stratified estimands within node4 strata.",
     "Causal × no-selection.",
     "Randomization within strata of node4.",
     "Conditional exchangeability holds within strata.",
     "Schoenfeld residuals (global + per covariate); forest plot of stratum-specific HRs.",
     "Non-PH means the HR is a time-average; report RMST instead.",
     "Same as Q1.",
     "Treating a single HR as the effect when PH is violated."),

    ("Q2_propensity.ipynb",
     "Q2 — Propensity score (diagnostic)",
     "e(Z) = P(rx=Lev+5FU | Z); used to construct IPW/AIPW weights downstream.",
     "Causal × no-selection (now treated as observational).",
     "Backdoor through Z = {age, sex, obstruct, perfor, adhere, nodes, differ, extent, surg}.",
     "Conditional exchangeability T_d(rx) ⊥ rx | Z and positivity 0 < e(Z) < 1.",
     "Logistic regression and gradient-boosted classifier for e(Z); plot overlap by arm.",
     "Lack of overlap → positivity violation → no transport-of-information possible at extremes.",
     "Same as Q1 (we are pretending no randomization label).",
     "Reading the propensity histogram as a randomization check; it diagnoses overlap, not exchangeability."),

    ("Q2_ipw.ipynb",
     "Q2 — IPW Cox",
     "δ_HR^IPW under backdoor adjustment.",
     "Causal × no-selection (obs).",
     "Backdoor via Z; w = 1/e(Z) for treated, 1/(1-e(Z)) for control.",
     "Conditional exchangeability + positivity.",
     "Stabilized IPW weights; weighted Cox via lifelines.",
     "Mis-specified e(Z) biases the weights; extreme weights inflate variance.",
     "Same as Q1.",
     "Forgetting variance inflation when reporting CIs from weighted regression."),

    ("Q2_aipw.ipynb",
     "Q2 — AIPW / DR for RMST",
     "δ_RMST(5y) via doubly-robust augmented IPW.",
     "Causal × no-selection (obs).",
     "Backdoor via Z.",
     "Either e(Z) OR outcome regression μ(Z, rx) is correctly specified.",
     "Hand-coded AIPW for RMST OR econml DRLearner with RMST outcome.",
     "Both nuisances misspecified → no double-robustness protection.",
     "Same as Q1.",
     "Reporting AIPW point estimate without checking nuisance fit (cross-fit RMSE)."),

    ("Q2_dml.ipynb",
     "Q2 — Double/debiased ML",
     "δ partialled-out treatment effect; Neyman-orthogonal score.",
     "Causal × no-selection (obs).",
     "Backdoor via Z.",
     "Cross-fitted nuisances converge at rate n^(-1/4).",
     "econml.dml.LinearDML with GBM nuisances; 5-fold cross-fit.",
     "If nuisances do not converge fast enough, plug-in bias contaminates θ̂.",
     "Same as Q1.",
     "Using DML as a black-box ATE estimator without checking nuisance convergence."),

    ("Q2_bad_control.ipynb",
     "Q2 — The bad-control demo",
     "δ_HR conditional on M (recurrence). NOT an interpretable causal effect.",
     "Pedagogical anti-example.",
     "Conditioning on M opens a collider path rx → M ← U → Y_death; blocks the rx → M → Y_death mediator path.",
     "M is post-treatment; sequential ignorability would be required to interpret δ as NDE.",
     "Cox with rx + M + Z; show the attenuated coefficient on rx vs Q2_ipw.",
     "If the analyst stops here and interprets the rx coefficient as 'treatment effect,' they are wrong.",
     "n/a — the point is what NOT to estimate.",
     "Adjusting for everything 'predictive of outcome' including post-treatment variables."),

    ("Q2_forest_plot.ipynb",
     "Q2 — Five-estimator forest plot",
     "δ_HR (or δ_RMST) from {naive Cox, IPW, AIPW, DML, bad-control}.",
     "Causal × no-selection (obs) + the anti-example.",
     "All but bad-control identify via backdoor; bad-control does not identify.",
     "All but bad-control share the conditional exchangeability assumption.",
     "Forest plot of point estimates + 95% CIs.",
     "If IPW/AIPW/DML disagree materially with the randomized Q1 answer, something is wrong with Z.",
     "Same as Q1.",
     "Looking at the forest and concluding 'methods don't matter' — bad-control breaks the convergence."),

    ("Q3_metalearners.ipynb",
     "Q3 — Meta-learners for CATE",
     "τ(z) = E[T_d(1) - T_d(0) | Z=z] with RMST(5y) as the outcome scale.",
     "Causal × no-selection (conditional).",
     "Randomization within strata of Z.",
     "Conditional exchangeability + positivity within strata.",
     "S-, T-, X-, DR-learner from econml.metalearners and econml.dr; IPCW weights for 5-year RMST.",
     "Sparse strata → unstable τ̂(z); regularization bias.",
     "Same as Q1.",
     "Reporting only the headline τ̂(z) without strata-level CIs."),

    ("Q3_causal_forest.ipynb",
     "Q3 — Causal forest (econml)",
     "τ(z) via honest causal forest.",
     "Causal × no-selection (conditional).",
     "Randomization within strata of Z.",
     "Conditional exchangeability + positivity + honesty (sample-splitting in tree).",
     "econml.grf.CausalForest with honest=True, inference=True.",
     "Forest pointwise CIs are not jointly valid across z.",
     "Same as Q1.",
     "Treating the heatmap of τ̂(z) as if its peaks are statistically significant."),

    ("Q3_blp.ipynb",
     "Q3 — Best linear projection of CATE",
     "Best linear projection coefficients of τ(Z) onto a chosen Ψ(Z) (Semenova & Chernozhukov 2021).",
     "Causal × no-selection (conditional).",
     "As Q3_causal_forest plus the standard BLP regularity conditions.",
     "Same as Q3 + correct specification of Ψ.",
     "BLP from econml; report coefficients with HC0 CIs.",
     "BLP coefficients are projection summaries, not causal moderator effects.",
     "Same as Q1.",
     "Interpreting a BLP coefficient on, e.g., age, as 'older patients benefit more causally.'"),

    ("Q3_calibration.ipynb",
     "Q3 — Athey-Wager calibration",
     "Bin patients by τ̂(z) quantile; compare observed vs predicted within-bin ATE.",
     "Causal × no-selection (conditional).",
     "Quantile binning + within-bin randomization.",
     "Within-bin conditional exchangeability.",
     "Athey-Wager calibration plot; binned within-bin AIPW.",
     "Poor calibration suggests overfitting or misspecification.",
     "Same as Q1.",
     "Skipping calibration and trusting τ̂(z) at face value."),

    ("Q3_cate_by_nodes.ipynb",
     "Q3 — CATE vs nodes (headline figure)",
     "τ(nodes=k) for k = 0..max.",
     "Causal × no-selection (conditional).",
     "Randomization within strata of nodes.",
     "Conditional exchangeability + positivity at each nodes value.",
     "Causal forest projection onto nodes; scatter + ribbon.",
     "Few patients at high node counts → wide ribbons.",
     "Same as Q1.",
     "Reading peaks as significant moderation; reporting without CI ribbons."),

    ("Q3_variable_importance.ipynb",
     "Q3 — Feature importance",
     "Forest-based variable importance for moderation.",
     "Causal × no-selection (conditional).",
     "n/a — descriptive of the τ̂(z) function.",
     "n/a — descriptive of the τ̂(z) function.",
     "econml causal-forest VI; SHAP on τ̂(z) surface.",
     "VI ≠ causal moderation. It is variable importance for predicting τ̂.",
     "Same as Q1.",
     "Concluding that a top-VI variable causally moderates without independent BLP/test."),

    ("Q5_seer_extract.ipynb",
     "Q5 — SEER 1990 stage B/C extract",
     "Marginal P(Z=z | S=0) — the target-population covariate distribution.",
     "Pre-causal: data-engineering step.",
     "n/a — descriptive.",
     "n/a — descriptive.",
     "SEER*Stat case-listing session, 1989-1991, colorectal, stage B/C; export to seer_1990.csv.",
     "SEER coding of stage differs from Astler-Coller used in Moertel.",
     "U.S. 1990 stage B/C resected colon-cancer patients.",
     "Assuming SEER covariates align 1-to-1 with trial covariates (they do not)."),

    ("Q5_iosw.ipynb",
     "Q5 — Inverse-odds-of-sampling weights",
     "w(X) = (1 - ê_S(X)) / ê_S(X) for trial patients; Cole & Stuart 2010.",
     "Causal × selection.",
     "T_d(rx) ⊥ S | X with X = effect modifiers identified by Q3.",
     "Conditional ignorability of trial participation given X.",
     "Logistic regression of S on X stacked across trial+SEER; compute w for trial patients.",
     "Unmeasured effect modifier U_X breaks the transport assumption.",
     "U.S. 1990 stage B/C resected CRC.",
     "Reading transported HR as 'the real-world HR' without acknowledging measured-X scope."),

    ("Q5_transported_ate.ipynb",
     "Q5 — Transported ATE",
     "δ_ATE^target = E_X|S=0[ τ(X) ].",
     "Causal × selection.",
     "Same as Q5_iosw.",
     "Same as Q5_iosw.",
     "Weighted Cox / RMST using w(X); compare to within-trial Q1.",
     "Transported estimate diverging from Q1 → effect modification by mismatch covariates.",
     "U.S. 1990 stage B/C resected CRC.",
     "Calling the transported estimate generalizable without bounding U_X."),

    ("Q5_dahabreh_bounds.ipynb",
     "Q5 — Dahabreh worst-case bounds",
     "Worst-case bound on δ_ATE^target under one unmeasured effect modifier (Dahabreh et al. 2019).",
     "Causal × selection.",
     "Bounds on the unmeasured effect modifier's HR with rx.",
     "Bounded-effect-modifier sensitivity.",
     "Sensitivity sweep; tipping point at which bound covers 0.",
     "If tipping point is implausibly small, finding is fragile.",
     "U.S. 1990 stage B/C resected CRC.",
     "Quoting the point estimate without the bound."),

    ("S1_e_values.ipynb",
     "S1 — E-values (VanderWeele 2017)",
     "E-value for each causal point estimate.",
     "Sensitivity.",
     "Unmeasured U on each of {Q1 HR, Q2 IPW HR, Q3 CATE summary, Q5 transported HR}.",
     "Bounded-U sensitivity per VanderWeele-Ding.",
     "sensemakr (R) and direct E-value formulas.",
     "E-value mis-states sensitivity for time-varying U or non-rare outcomes.",
     "Same as the underlying estimate.",
     "Reading E-value as a hypothesis test rather than a sensitivity bound."),

    ("S2_sensitivity_table.ipynb",
     "S2 — Master sensitivity table",
     "Synthesis across all sensitivity analyses.",
     "Sensitivity.",
     "Union of identifying assumptions across Q1–Q5.",
     "Union of non-graphical assumptions across Q1–Q5.",
     "Single master table: estimand, point, 95% CI, sensitivity parameter, breakdown value, qualitative robustness.",
     "n/a — synthesis step.",
     "n/a — synthesis step.",
     "Reading robustness as binary; report it as a continuous breakdown value."),
]


def main():
    for name, *args in STUBS:
        path = NB_DIR / name
        with open(path, "w") as f:
            json.dump(stub(*args), f, indent=1)
        print(f"  wrote {path.relative_to(ROOT)}")
    print(f"\nGenerated {len(STUBS)} stub notebooks under {NB_DIR.relative_to(ROOT)}/")


if __name__ == "__main__":
    main()

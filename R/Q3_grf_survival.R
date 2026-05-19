# =============================================================================
# Q3 — Causal survival forest (gold-standard, R grf package)
# =============================================================================
# Estimand:       tau(z) = E[T_d(rx=Lev+5FU) - T_d(rx=Obs) | Z=z], RMST(5y)
# Cell:           Causal x no-selection (conditional)
# Identifies via: Randomization within strata of Z (graphical), conditional
#                 exchangeability + positivity (non-graphical).
# Estimator:      grf::causal_survival_forest with horizon = 5 * 365.25 days.
# Failure mode:   Sparse strata at high node counts; pointwise CIs only.
# Population:     Stage B/C resected colon cancer patients (Moertel trial).
# Naive mistake:  Treating peaks of tau(z) as significant moderation.
#
# Compared to econml.grf.CausalForest, grf's causal_survival_forest natively
# handles censoring on the time scale, which is the cleanest specification for
# RMST-style estimands.
# =============================================================================

# Filled in Week 4.

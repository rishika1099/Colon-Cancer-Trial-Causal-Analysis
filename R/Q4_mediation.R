# =============================================================================
# Q4 — Mediation analysis (Imai-Keele-Tingley natural effects)
# =============================================================================
# Estimand:
#   NIE = E[T_d(rx=1, M(1)) - T_d(rx=1, M(0))]
#   NDE = E[T_d(rx=1, M(0)) - T_d(rx=0, M(0))]
#   TE  = NDE + NIE
# Cell:             Causal x no-selection (mediator decomposition).
# Identifies via:
#   (i)  Treatment ignorability given baseline Z (randomization, satisfied).
#   (ii) Mediator ignorability given treatment + baseline Z (the assumption
#        that bites — unobserved frailty U could affect both M and Y).
# Estimator:        mediation::mediate(), bootstrap=1000, with
#                   outcome model:  coxph(Surv(t_death, status) ~ rx + M + Z)
#                   mediator model: glm(M ~ rx + Z, family=binomial)
# Failure mode:     U->M and U->Y break sequential ignorability; bound via
#                   mediation::medsens (Imai-rho).
# Population:       Stage B/C resected colon cancer patients (Moertel trial).
# Naive mistake:    Calling the M-adjusted Cox coefficient on rx the "direct
#                   effect." That is NOT NDE — see Q2_bad_control.
# =============================================================================

# Filled in Week 5–6.

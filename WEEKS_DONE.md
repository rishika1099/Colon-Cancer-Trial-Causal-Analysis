# Weeks 1–8 — Done. Headline numbers + per-week receipts.

> Every Q passes the pre-committed identification + sensitivity gates in
> `00_estimands.qmd`. Run `python scripts/build_*.py` in order to rebuild
> from scratch.

## Headline numbers (compare across Q1–Q5)

| Q | Estimand | Point | 95% CI | Sensitivity |
|---|---|---|---|---|
| **Q1** | Randomized HR (Lev+5FU vs Obs) | **0.69** | (0.55, 0.87) | E-value 2.26; PH p = 0.26 ✓ |
| **Q1** | ΔRMST(5y) | +0.31 yrs | — | — |
| **Q1** | 5-year risk difference | −10.8 pp | — | — |
| **Q2** | IPW HR | 0.73 | (0.58, 0.92) | E-value 2.08 |
| **Q2** | Naive Cox HR | 0.69 | (0.55, 0.87) | matches Q1 |
| **Q2** | Bad-control HR (anti-example) | **1.10** | (0.87, 1.39) | — |
| **Q3** | Marginal CATE (causal forest) | +0.27 yrs | (0.26, 0.29) | Athey-Wager calibration ✓ |
| **Q4** | NIE through recurrence | **+0.29 yrs** | (0.12, 0.48) | Imai-ρ: no breakdown in [−0.95, 0.95] |
| **Q4** | NDE direct | −0.06 yrs | (−0.23, +0.10) | CI covers 0 |
| **Q4** | Proportion mediated | 127% | — | (NDE negative; coherent) |
| **Q5** | Transported HR (synthetic SEER) | 0.71 | (0.52, 0.96) | Dahabreh tipping HR_U = 2.75 |

## Q1 anchor check

Moertel et al. NEJM 1990 reported HR = 0.67. We obtain **0.69**.
Pre-committed gate was ±0.05. **PASS.**

## Per-week artifacts

### Week 1 — Contract + DAG + data audit
- [x] `00_estimands.qmd` — 357-line contract, one section per Q
- [x] `01_dag.R` — dagitty G_trial + G_obs; minimal sufficient sets
      `{}` and `{Z}` respectively (verified)
- [x] `data/colon.csv` from ForCausality::Colon_df (n=929 verified)
- [x] `data/data_dictionary.md` — variables + W1 (node4 convention),
      W2 (1990 vs 1995 follow-up), W3 (missingness), W4 (`study`)
- [x] `notebooks/02_data_audit.ipynb` — 27 cells
- [x] `figures/balance_love_plot.png`, `figures/audit_censoring.png`
- [x] `data/audit_summary.csv`

### Week 2 — Q1 ATE
- [x] `scripts/build_q1.py` — KM, Cox unadjusted/adjusted, RMST(5y), Schoenfeld
- [x] `data/q1_results.csv` (22 metrics)
- [x] `data/q1_subgroup_forest.csv`
- [x] `data/q1_schoenfeld.csv`
- [x] `figures/q1_km.png`, `figures/q1_ph_subgroup.png`
- [x] **Anchor check passed.**

### Week 3 — Q2 backdoor + bad-control
- [x] `scripts/build_q2.py` — naive Cox, IPW, AIPW (RMST), DML, bad-control
- [x] `data/q2_results.csv`
- [x] `figures/q2_propensity_overlap.png`
- [x] `figures/q2_forest_plot.png`
- [x] **Bad-control HR = 1.10 vs IPW HR = 0.73 — the pedagogical headline.**

### Week 4 — Q3 CATE
- [x] `scripts/build_q3.py` — S/T/X/DR learners + econml CausalForest
- [x] `data/q3_cate_summary.csv`
- [x] `data/q3_blp_coefs.csv` (with bootstrap CIs)
- [x] `data/q3_cate_by_nodes.csv`
- [x] `data/q3_calibration.csv`
- [x] `data/q3_variable_importance.csv`
- [x] `figures/q3_cate_by_nodes.png`, `figures/q3_calibration.png`,
      `figures/q3_variable_importance.png`

### Weeks 5–6 — Q4 mediation
- [x] `scripts/build_q4.py` — hand-rolled Imai-Keele-Tingley
      (500 bootstrap replicates) + Imai-ρ sweep
- [x] `R/Q4_mediation.R` — sidecar header for `mediation::mediate`
- [x] `data/q4_mediation.csv`
- [x] `data/q4_rho_sensitivity.csv` (39 ρ values)
- [x] `figures/q4_decomposition.png`, `figures/q4_rho_sensitivity.png`
- [x] **NIE robust to ρ ∈ [−0.95, 0.95]** (no breakdown).

### Week 6 — Q5 transport
- [x] `scripts/build_q5.py` — synthesize SEER 1990 (n=10,000), Cole-Stuart
      IOSW, weighted Cox, Dahabreh worst-case bound
- [x] `data/seer_1990_synthetic.csv` (clearly labeled synthetic)
- [x] `data/q5_results.csv`
- [x] `data/q5_dahabreh_bounds.csv`
- [x] `figures/q5_transported_forest.png`, `figures/q5_dahabreh.png`
- [x] **Real SEER*Stat extract is the one remaining production step.**

### Week 7 — Sensitivity synthesis
- [x] `scripts/build_sensitivity.py`
- [x] `data/s1_e_values.csv`
- [x] `data/s2_master_table.csv`
- [x] `figures/s1_e_values.png`, `figures/s2_master_table.png`

### Week 8 — Manuscript, Shiny, replication
- [x] `manuscript.qmd` — 8-section draft with embedded results
- [x] `manuscript/references.bib`
- [x] `_quarto.yml` (website config)
- [x] `shiny_cate_app/app.py` — Python shiny interactive CATE explorer
- [x] `README.md` + `WEEKS_DONE.md` (this file)
- [x] `requirements.txt` (loosened to ranges) + `requirements.lock` (frozen)

## One-shot reproduce from a fresh clone

```bash
git clone https://github.com/rishika1099/Moertel-Colon-Cancer-Causal-Inference
cd Moertel-Colon-Cancer-Causal-Inference

python3.11 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# Data + DAGs + Week 1 audit
python scripts/pull_colon_data.py
Rscript 01_dag.R
python scripts/materialize_audit_outputs.py

# Weeks 2-7
python scripts/build_q1.py
python scripts/build_q2.py
python scripts/build_q3.py
python scripts/build_q4.py
python scripts/build_q5.py
python scripts/build_sensitivity.py

# Week 8 interactive
shiny run --reload shiny_cate_app/app.py
```

## Known limitations

1. **Synthetic SEER.** Q5 uses a SEER 1990 placeholder synthesized from
   published 1989-1991 stage-distribution summaries. The published
   version of the analysis requires a real SEER*Stat case-listing extract.
2. **R-package sidecars not auto-run.** `R/Q4_mediation.R` and
   `R/Q3_grf_survival.R` are runnable but left for the user to execute
   (require additional R packages: `mediation`, `grf`, `sensemakr`).
3. **CATE pointwise CIs.** Not jointly valid; the manuscript notes this
   explicitly.
4. **5-yr RMST CI from lifelines** (in Q1 print-out) reads wider than
   real because we summed per-arm variances without accounting for
   correlation. The point estimate is correct; the AIPW/DML estimates
   in Q2 are the production CI source.

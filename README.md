# Moertel Colon Cancer Causal Inference

Re-analysis of the Moertel et al. NEJM 1990 trial of adjuvant levamisole +
5-fluorouracil in resected stage B/C colon cancer (`n = 929`) using
five causal-inference methods.

The dataset is the well-known `colon` dataset from the R `survival` package
(also redistributed as `ForCausality::Colon_df`). Each patient has two
rows: one for recurrence (`etype = 1`) and one for death (`etype = 2`).
The primary analytic frame is the death-row frame; the recurrence row is
the mediator in Q4.

The trial's published hazard ratio is approximately `0.67`. The
re-analysis recovers `0.69` (95% CI `0.55 – 0.87`) and then layers four
additional structural questions on top of it.

## Results

| Q | Estimand | Estimate | 95% CI |
|---|---|---|---|
| Q1 | Randomized Cox HR (Lev+5FU vs Obs) | 0.690 | (0.55, 0.87) |
| Q1 | 5-year RMST contrast | +0.31 yrs | — |
| Q1 | 5-year risk difference | −10.8 pp | — |
| Q2 | IPW Cox HR | 0.730 | (0.58, 0.92) |
| Q2 | AIPW 5-yr RMST contrast | +0.281 yrs | (+0.09, +0.48) |
| Q2 | DML 5-yr RMST contrast | +0.179 yrs | (−0.06, +0.42) |
| Q2 | Bad-control Cox HR (conditions on M) | 1.103 | (0.87, 1.39) |
| Q3 | Marginal CATE (causal forest, RMST) | +0.27 yrs | (+0.26, +0.29) |
| Q4 | Natural indirect effect (NIE) | +0.292 yrs | (+0.12, +0.48) |
| Q4 | Natural direct effect (NDE) | −0.062 yrs | (−0.23, +0.10) |
| Q5 | Transported HR (synthetic SEER 1990) | 0.707 | (0.52, 0.96) |

Sanity checks against the published trial: 5-year survival 52.6 / 53.5 /
63.4% for Obs / Lev / Lev+5FU (Moertel reported about 55% / 55% / 65%);
levamisole-alone Cox HR `0.974` with CI crossing 1 (the original paper's
secondary null result); NDE + NIE = TE to numerical precision.

The Q2 bad-control row is not a competitor estimator. It is included to
show what happens when a post-treatment mediator (recurrence) is added
to a Cox model that is otherwise correctly specified: the apparent HR
moves from `0.69` to `1.10`. This is the size of the error a careless
adjustment can produce on real survival data.

## The five questions

| Q | Question | Identification |
|---|---|---|
| Q1 | Randomized ATE | Randomization severs `Z → rx` |
| Q2 | Back-door identification (forget randomization) | Adjust on `Z` in `G_obs` |
| Q3 | Heterogeneous effects, `τ(z)` | Randomization within strata of `Z` |
| Q4 | Mediation through recurrence | Sequential ignorability (Imai-Keele-Yamamoto) |
| Q5 | Transport to SEER 1990 stage B/C | Cole–Stuart inverse-odds-of-sampling |

Two DAGs are committed and verified with `dagitty`:

- `G_trial`: randomization severs `Z → rx`. Minimal sufficient adjustment
  set for `rx → Y_death` is the empty set.
- `G_obs`: counterfactual world where `Z → rx` is restored. Minimal
  sufficient adjustment set is `{Z}`.

<p align="center">
  <img src="figures/dag_trial.png" width="48%" alt="G_trial DAG"/>
  <img src="figures/dag_observational.png" width="48%" alt="G_observational DAG"/>
</p>

## Report

The IEEE conference-format write-up of the project is in
[`report/ieee_report.pdf`](report/ieee_report.pdf) (5 pages, US Letter).
Source: [`report/ieee_report.tex`](report/ieee_report.tex). Build it
locally with `cd report && tectonic ieee_report.tex` (or `pdflatex` if
you have TeX Live).

## Repository layout

```
.
├── 00_estimands.qmd              estimand sheet (one section per Q)
├── 01_dag.R                      dagitty DAGs, SVG/PNG export
├── manuscript.qmd                long-form manuscript draft
├── manuscript/references.bib     bibliography
├── report/                       IEEE conference paper (.tex + .pdf)
├── index.qmd                     Quarto site landing page
├── _quarto.yml                   Quarto website config
│
├── data/
│   ├── colon.csv                 929 patients x 16 cols, from ForCausality
│   ├── seer_1990_synthetic.csv   synthesized SEER target for Q5
│   ├── data_dictionary.md        per-variable docs, including data warts
│   └── q{1..5}_*.csv, s{1,2}_*   results tables produced by build scripts
│
├── scripts/
│   ├── pull_colon_data.py        rpy2 -> ForCausality::Colon_df -> data/colon.csv
│   ├── build_q1.py ... build_q5.py one script per question
│   ├── build_sensitivity.py      E-values + master sensitivity table
│   ├── build_audit_notebook.py   programmatic builder for the audit notebook
│   ├── materialize_audit_outputs.py
│   └── make_stub_notebooks.py
│
├── notebooks/
│   ├── 02_data_audit.ipynb       executed balance + missingness + warts audit
│   └── Q1_*, Q2_*, Q3_*, Q5_*, S{1,2}_*   estimand-headed stubs
│
├── R/
│   ├── Q3_grf_survival.R         grf causal_survival_forest sidecar
│   └── Q4_mediation.R            mediation::mediate + medsens sidecar
│
├── shiny_cate_app/app.py         interactive CATE explorer (Python shiny)
├── figures/                      all manuscript and report figures
├── requirements.txt              human-readable pins (ranges)
├── requirements.lock             pip freeze for byte-exact reproducibility
└── LICENSE                       MIT
```

## Reproducibility

```bash
git clone https://github.com/rishika1099/Moertel-Colon-Cancer-Causal-Inference
cd Moertel-Colon-Cancer-Causal-Inference

# Python 3.11 environment
python3.11 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# R (only required for the DAG drawings and the optional mediation/grf sidecars)
Rscript -e 'install.packages(c("dagitty","svglite","ForCausality","survival"))'

# Build pipeline
python scripts/pull_colon_data.py            # data
Rscript 01_dag.R                              # DAGs
python scripts/materialize_audit_outputs.py   # audit artifacts
python scripts/build_q1.py
python scripts/build_q2.py
python scripts/build_q3.py
python scripts/build_q4.py
python scripts/build_q5.py
python scripts/build_sensitivity.py

# Interactive CATE explorer
shiny run --reload shiny_cate_app/app.py

# Optional Quarto site / PDF
quarto render
```

Random seeds are fixed in every script (`np.random.seed(42)` in Python,
`set.seed(42)` in R).

## Methods summary

**Q1 — randomized ATE.** Kaplan-Meier with log-rank, Cox proportional
hazards (unadjusted and covariate-adjusted), 5-year restricted mean
survival time, Schoenfeld global PH test, subgroup forest by lymph-node
burden.

**Q2 — back-door identification.** Logistic and gradient-boosted
propensity scores with overlap diagnostic; stabilized IPW Cox;
hand-coded AIPW for 5-year RMST with IPCW censoring weights;
`econml.LinearDML` with 5-fold cross-fitting; regression-adjusted Cox;
and a deliberate bad-control Cox conditioning on the recurrence
mediator.

**Q3 — heterogeneous treatment effects.** S-, T-, X-, and DR-learners
(via `econml.metalearners` and `econml.dr`), and an honest causal forest
(`econml.grf.CausalForest`, `honest=True`, `inference=True`).
Best-linear projection of `τ̂(Z)` onto `Z` with 200-replicate bootstrap
CIs. Athey-Wager binned calibration at five quantile bins.

**Q4 — mediation.** Python implementation of the Imai-Keele-Tingley
(2010) Algorithm 1: logistic mediator model, gradient-boosted outcome
model on IPCW-weighted 5-year RMST, 500 bootstrap replicates. Imai-`ρ`
sensitivity sweep over 39 values in `[-0.95, 0.95]`. R
`mediation::mediate` is provided as a sidecar for cross-validation.

**Q5 — transportability.** Synthetic SEER 1989–1991 stage B/C target
(`n = 10,000`) generated from published stage-distribution summaries.
Cole-Stuart inverse-odds-of-sampling weights with a logistic
propensity-of-trial model; weighted Cox on the trial sample. Dahabreh et
al. (2019) worst-case bound sweep over an unmeasured effect modifier.

**Sensitivity synthesis.** VanderWeele-Ding E-values for every HR
estimate, plus a master sensitivity table mapping each estimand to its
point estimate, CI, sensitivity parameter, breakdown value, and a
qualitative robustness label.

## Notes on scope

- The Q5 SEER target is synthetic. Producing it from a real SEER\*Stat
  case-listing extract is a one-CSV change.
- The R sidecars (`R/Q4_mediation.R`, `R/Q3_grf_survival.R`) are
  runnable but require additional R packages (`mediation`, `grf`,
  `sensemakr`). The Python implementations are the source of the numbers
  reported in this README and the IEEE write-up.
- The Q3 confidence intervals are pointwise. They are not jointly valid
  across the covariate surface.
- The Q4 proportion mediated exceeds 100% because the NDE is negative
  while the NIE is positive. NDE + NIE = TE holds exactly.

## Citation

Original trial:

> Moertel CG, Fleming TR, Macdonald JS, *et al.* Levamisole and
> fluorouracil for adjuvant therapy of resected colon carcinoma.
> *N Engl J Med* 1990;322(6):352–358.
> [doi:10.1056/NEJM199002083220602](https://doi.org/10.1056/NEJM199002083220602)

This re-analysis:

```bibtex
@misc{mamidibathula2026moertel,
  title  = {A Causal-Inference Re-analysis of the Moertel 1990 Adjuvant Colon Cancer Trial},
  author = {Mamidibathula, Rishika},
  year   = {2026},
  howpublished = {\url{https://github.com/rishika1099/Moertel-Colon-Cancer-Causal-Inference}}
}
```

## License

MIT. See [`LICENSE`](LICENSE).

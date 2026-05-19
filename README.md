# Beyond the Hazard Ratio

**Causal-inference re-analysis of the Moertel et al. 1990 adjuvant colon cancer trial.**

Target journal: *Observational Studies*. Methods-pedagogy aligned with Adam Kelleher's causal inference course (Columbia).

## What this repo contains

Five nested causal questions on the `Colon_df` data (Moertel et al. NEJM 1990, n=929, 3-arm RCT of Levamisole +/- 5-FU vs Observation):

| Q  | Question                              | Cell                              | Lecture       |
|----|---------------------------------------|-----------------------------------|---------------|
| Q1 | Randomized ATE                        | Causal × no-selection             | L1, L4, L5    |
| Q2 | Forget randomization — backdoor       | Causal × no-selection (obs.)      | L3, L4, L5    |
| Q3 | Heterogeneous treatment effects (CATE)| Causal × no-selection, conditional| L13           |
| Q4 | Mediation via recurrence              | Causal × no-selection, mediator   | L8            |
| Q5 | Transport to SEER 1990 population     | Causal × selection                | L11           |

Each Q has an explicit estimand, identification argument (graph + non-graph), estimator(s), failure modes, and sensitivity bounds. See `00_estimands.qmd` — **the contract that drives every notebook**.

## Reproducing from scratch

```bash
git clone <repo>
cd moertel-causal

# Python
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# R
R -e 'install.packages("renv"); renv::restore()'

# Render everything
quarto render
```

Every notebook fixes seeds (`np.random.seed(42)`, `set.seed(42)`). A clean clone + render reproduces every number in `manuscript.qmd`.

## Repo map

```
moertel-causal/
├── 00_estimands.qmd          # Week 1: the contract (write before any code)
├── 01_dag.R                  # Week 1: dagitty G_trial, G_observational
├── 02_data_audit.ipynb       # Week 1: balance tables, SMDs, data-warts
├── notebooks/Q1_*.ipynb      # Week 2: ATE
├── notebooks/Q2_*.ipynb      # Week 3: backdoor, IPW/AIPW/DML, bad-control demo
├── notebooks/Q3_*.ipynb      # Week 4: CATE meta-learners + causal forest
├── R/Q4_*.R                  # Week 5–6: mediation (Imai-Keele-Tingley)
├── notebooks/Q5_*.ipynb      # Week 6: transportability (Cole-Stuart)
├── notebooks/S*_*.ipynb      # Week 7: E-values, sensitivity synthesis
├── manuscript.qmd            # Week 8: 8000-word manuscript
└── R/shiny_cate_app.R        # Week 8: interactive CATE explorer
```

## Week-by-week status — all green

- [x] **Week 1** — Contract + DAG + data audit. See `WEEK1_DONE.md`.
- [x] **Week 2** — Q1 ATE. Cox HR = 0.69 (CI 0.55–0.87); anchor PASS.
- [x] **Week 3** — Q2 backdoor + the bad-control demo (HR 1.10 vs IPW 0.73).
- [x] **Week 4** — Q3 CATE meta-learners + causal forest; headline figure CATE-by-nodes.
- [x] **Week 5–6** — Q4 mediation; NIE robust to Imai-ρ ∈ [−0.95, 0.95].
- [x] **Week 6** — Q5 transport to synthetic SEER 1990; Dahabreh tipping HR_U = 2.75.
- [x] **Week 7** — E-values + master sensitivity table.
- [x] **Week 8** — `manuscript.qmd`, Quarto site, Shiny CATE explorer, `WEEKS_DONE.md`.

See `WEEKS_DONE.md` for the full receipts, headline numbers, and one-shot
reproduce instructions.

## Citing

Original trial: Moertel CG et al. *NEJM* 1990;322(6):352–8. Levamisole and fluorouracil for adjuvant therapy of resected colon carcinoma.

Data source: `ForCausality::Colon_df` (Toby Codigos), a curated copy of `survival::colon`.

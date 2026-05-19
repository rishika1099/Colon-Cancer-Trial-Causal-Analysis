# Week 1 — Done. Review before greenlighting Week 2.

> The estimand document is the contract. Read it carefully before signing
> off. If anything in `00_estimands.qmd` is wrong, every downstream notebook
> inherits the error.

## What got built

### Repo scaffolding

- [x] Directory layout per the spec (`R/`, `notebooks/`, `data/`,
      `figures/`, `manuscript/`, `scripts/`).
- [x] `.gitignore`, `requirements.txt` (pinned), `_quarto.yml`,
      `index.qmd`, `README.md`.
- [x] **21 stub notebooks** in `notebooks/`, one per Q-step. Every stub
      opens with the estimand-first 7-field header (estimand, cell,
      identification graph, identification non-graph, estimator, failure
      mode, target population, naive-reader mistake) so Week 2+ can fill
      in analysis cells without rewriting the rigor scaffolding.
- [x] R-script stubs `R/Q3_grf_survival.R`, `R/Q4_mediation.R` with
      header-level estimand commitments.

### The contract

- [x] **`00_estimands.qmd`** — full draft, one section per Q (Q1–Q5),
      plus preliminaries, the 2×2 cell map, the two DAGs verbalized,
      sensitivity synthesis, and a closing "what this document does *not*
      commit to." Uses Adam's Kelleher-course vocabulary throughout
      (Pearl's ladder, $\delta_{\text{naive}} \to \delta$ under
      $Y_d^{(rx)} \perp rx$, back-door criterion, sequential
      ignorability, transportability per Cole-Stuart).
- [x] **`01_dag.R`** — `dagitty` for $G_{\text{trial}}$ and
      $G_{\text{observational}}$. Verified by dagitty:
      - $G_{\text{trial}}$ minimal adjustment set: $\emptyset$ ✓
        (randomization)
      - $G_{\text{observational}}$ minimal adjustment set: $\{Z\}$ ✓
      - Implied independencies: $S \perp rx$, $Z \perp rx$,
        $M \perp S \mid Z$, $S \perp Y_{\text{death}} \mid Z$
- [x] SVG + PNG figures: `figures/dag_trial.{svg,png}`,
      `figures/dag_observational.{svg,png}`, and
      `figures/dag_adjustment_sets.txt` (the dagitty receipt).

### Data + audit

- [x] **`scripts/pull_colon_data.py`** — pulls
      `ForCausality::Colon_df` via rpy2 (with an Rscript-subprocess
      fallback for fresh clones). Verifies 1,858 rows, 929 unique patient
      IDs, 929 etype==2 rows, `rx` levels = {Obs, Lev, Lev+5FU}. Writes
      `data/colon.csv`.
- [x] **`notebooks/02_data_audit.ipynb`** — full Week 1 audit, 27
      cells, with the estimand-first header (descriptive cell) plus:
      1. Load & shape
      2. Arm sizes
      3. SMD balance table + Love plot
      4. `nodes` vs `node4` forensic check (W1)
      5. Missingness (W3)
      6. Censoring + event summary
      7. Crude 5-year death rate by arm (Q1 anchor)
      8. Summary write-out
      9. What this audit does *not* establish
- [x] **`data/data_dictionary.md`** — per-variable descriptions plus a
      "Data warts" section documenting W1 (`node4` threshold), W2 (1990
      vs 1995 follow-up), W3 (missingness), W4 (vestigial `study` column).
- [x] **`data/audit_summary.csv`** — compact metrics table.
- [x] **`figures/balance_love_plot.png`** + **`audit_censoring.png`** —
      figures referenced by the manuscript draft.

## Numbers to verify before greenlighting Week 2

These are the values the audit and DAG scripts produced. **If any of these
look wrong, stop and investigate.**

| Quantity | Value |
|---|---|
| Raw rows (both etype) | 1,858 |
| Unique patients | 929 |
| Patients in death-row frame | 929 |
| Arm sizes (Obs / Lev / Lev+5FU) | 315 / 310 / 304 |
| Overall death rate | 48.7% |
| Overall recurrence rate | 50.4% |
| Max |SMD| across arms | 0.126 (< 0.15 — well balanced) |
| Missing `nodes` | 18 (1.9%) |
| Missing `differ` | 23 (2.5%) |
| `nodes`/`node4` mismatches under `>=` | 90 |
| `nodes`/`node4` mismatches under `>` (true errors) | **12** |
| 5-year death rate — Obs | 47.3% |
| 5-year death rate — Lev | 46.5% |
| 5-year death rate — Lev+5FU | 36.5% |
| **5-year crude risk difference (Lev+5FU vs Obs)** | **−10.8 percentage-points** |

The −10.8pp 5-year risk difference is the **Q1 anchor**. The Moertel 1990
HR of 0.67 is consistent with this under a Weibull-ish baseline hazard;
Q1's properly identified estimate should reproduce both.

## Adam-style rigor checklist for the estimand document

Walk through each section of `00_estimands.qmd` and confirm:

- [ ] **Q1.** Estimand is written in $do(\cdot)$ or
      potential-outcomes notation, *not* as a regression coefficient.
      Cell labeled. Both graph and non-graph identification arguments
      stated. PH failure mode acknowledged (Schoenfeld + RMST companion).
      The 0.67 ground-truth anchor is in the section.
- [ ] **Q2.** Backdoor argument cites `01_dag.R`'s adjustment set $\{Z\}$.
      The five-estimator forest (incl. bad-control demo) is committed
      to. Sequential ignorability vs the bad-control coefficient is *not*
      conflated — the difference is explicit.
- [ ] **Q3.** CATE estimand is the conditional expectation, not the
      individual treatment effect. The pointwise-vs-joint CI caveat
      appears. Variable importance is *not* called "moderation."
- [ ] **Q4.** Imai-Keele-Yamamoto's *two* assumptions stated
      separately. The natural-effects formulation is distinguished from
      Adam's Lecture 8 front-door (front-door does *not* apply here
      because the direct `rx → Y_death` path exists). $\rho^*$
      breakdown thresholds (0.30 / 0.15) are committed to before seeing
      the data.
- [ ] **Q5.** Cole-Stuart transportability assumption is written as
      $\{Y_d^{(rx=r)}\}_r \perp S \mid X$, *not* as "weighted regression."
      The SEER-covariate-mismatch threat is named. Dahabreh bound is
      committed to.
- [ ] **Causal language sweep.** No claim of the form "X *causes* Y"
      appears without an identification argument upstream. No
      "associated with" language smuggled in for what should be a causal
      claim.
- [ ] **Scope language.** Every Q ends with an explicit target
      population statement; none says "the result generalizes."
- [ ] **Sensitivity commitments.** Every causal estimate has a
      pre-specified sensitivity bound listed in `§S1` / `§S2`.

## How to reproduce Week 1 from a fresh clone

```bash
cd moertel-causal

# 1. Data — uses Rscript subprocess (no rpy2 needed yet).
python3 scripts/pull_colon_data.py

# 2. DAGs.
Rscript 01_dag.R

# 3. Audit artifacts (figures + summary CSV) without spinning up Jupyter.
python3 scripts/materialize_audit_outputs.py

# 4. (Optional, with Jupyter env) full notebook with markdown narrative.
jupyter nbconvert --to notebook --execute --inplace notebooks/02_data_audit.ipynb
```

## What is *not* in Week 1 (deliberately deferred)

- No Kaplan-Meier curves, no Cox models, no HR estimates. Those are Q1.
- No propensity scores, no IPW. Those are Q2.
- No CATE work. That's Q3.
- No mediation. That's Q4.
- No SEER pull. That's Q5.
- No `manuscript.qmd` draft. That's Week 8.

## Open questions to resolve before Week 2

1. **PH check tolerance.** I committed `|Δ| ≤ 0.05` between our Q1 HR and
   Moertel's 0.67 as the debugging gate. Confirm or relax.
2. **SEER access.** Q5 requires SEER*Stat. Confirm you have an account /
   the case-listing extract is downloadable by Week 6.
3. **`mediation` R package or Python hand-roll first?** I committed to
   *both* (R via rpy2 *and* hand-rolled in Python for reviewer
   transparency). Confirm bandwidth — the hand-roll is ~100 lines.
4. **Lev-alone secondary contrasts.** I pre-specified that Lev-alone is
   reported but is not the headline. Confirm.
5. **Imai-ρ breakdown thresholds.** I committed 0.30 / 0.15 as the
   moderate / fragile cutoffs *before* seeing Q4 results. Confirm or
   tighten.

---

**Recommended next step.** Send `00_estimands.qmd` to your Adam-style
reviewer (or paste the Q1–Q5 sections into a chat). The reviewer's job is
to push back on language drift, hand-waving, and assumption omissions
*before* a single Cox model is fit. The contract is upstream; rigor here
prevents downstream rework.

When that review comes back and the document is locked, Week 2 (Q1) can
begin.

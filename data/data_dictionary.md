# Data dictionary — `data/colon.csv`

**Source.** `ForCausality::Colon_df` (CRAN), which is a verbatim copy of
`survival::colon` (R `survival` package). The underlying trial is Moertel
et al. NEJM 1990 (10.1056/NEJM199002083220602), with extended follow-up
through ~1995 (Moertel et al. Ann Intern Med 1995).

**Shape.** 1,858 rows × 16 columns. Two rows per patient — one for
recurrence (`etype == 1`) and one for death (`etype == 2`) — so the unique
patient count is 929.

**Primary analytic frames.**

| Frame | Rows | Used for |
|-------|------|----------|
| `df[df.etype == 2]` (death-row) | 929 | Q1, Q2, Q3, Q5 (primary outcome) |
| `df[df.etype == 1]` (recurrence-row) | 929 | Q4 (`M` = recurrence-as-mediator), Q2 bad-control demo |

---

## Variables

### Identifiers

| Variable | Type | Description |
|----------|------|-------------|
| `id`     | int  | Patient ID, 1–929. Stable across the two etype rows. |
| `study`  | int  | Always `1`. Legacy field from the multi-trial origin of the dataset. |
| `etype`  | int  | `1` = recurrence row; `2` = death row. Determines what `time` and `status` mean for that row. |

### Treatment (causal state)

| Variable | Type | Values | Description |
|----------|------|--------|-------------|
| `rx`     | str  | `Obs`, `Lev`, `Lev+5FU` | Randomized treatment arm. `Obs` = observation only (concurrent control). `Lev` = levamisole. `Lev+5FU` = levamisole + 5-fluorouracil. Headline contrast: `Lev+5FU` vs `Obs`. |

Arm sizes (death-row frame): Obs 315, Lev 310, Lev+5FU 304.

### Baseline covariates $Z$

| Variable   | Type    | Description | Notes |
|------------|---------|-------------|-------|
| `age`      | int     | Age at randomization, years (range 18–85). | Reported in years. |
| `sex`      | binary  | `1` = male, `0` = female. | ~52% male overall. |
| `obstruct` | binary  | `1` if tumor obstructed the colon at presentation. | ~19% prevalence. |
| `perfor`   | binary  | `1` if tumor perforated the bowel wall. | ~3% prevalence (rare). |
| `adhere`   | binary  | `1` if tumor adhered to adjacent organs. | ~14% prevalence. |
| `nodes`    | int     | Number of histologically positive lymph nodes detected at surgery. | **18 missing** in the death-row frame. Range 0–33. |
| `differ`   | int     | Tumor differentiation grade: `1` = well, `2` = moderate, `3` = poor. | **23 missing**. |
| `extent`   | int     | Local extent of tumor: `1` = submucosa, `2` = muscle, `3` = serosa, `4` = contiguous structures. | No missing. |
| `surg`     | binary  | `1` if time from surgery to randomization was *long* (>30 days). | ~27% prevalence. |
| `node4`    | binary  | **Encoded as `I(nodes > 4)`, not `I(nodes >= 4)`.** See "Data warts" below. |

### Outcomes

The same two column names (`time`, `status`) describe different things
depending on `etype`:

| Variable | When `etype==1` | When `etype==2` |
|----------|-----------------|-----------------|
| `time`   | Days from randomization to first recurrence or last recurrence-free follow-up. | Days from randomization to death or last live follow-up. |
| `status` | `1` if recurrence observed, `0` if censored. | `1` if death observed, `0` if censored. |

Recurrence rate (overall): 468 / 929 ≈ 50.4%.
Death rate (overall, full follow-up): 452 / 929 ≈ 48.7%.
Max follow-up: 3,329 days ≈ 9.1 years.

### Auxiliary

| Variable | Type | Description |
|----------|------|-------------|
| `S`      | added downstream | Trial-enrollment indicator. Equals `1` for every Moertel row; will equal `0` for SEER rows when those are stacked in Q5. Not in the raw CSV. |

---

## Data warts (documented per project rigor rule #8)

### W1. `node4` vs `nodes` — threshold convention

The variable name suggests `node4 = I(nodes >= 4)`, but the dataset encodes
`node4 = I(nodes > 4)` (i.e., strictly more than 4 positive nodes,
equivalently 5 or more). The 90 apparent disagreements vanish to 12 under
the correct convention.

**Twelve genuine inconsistencies remain** even under `nodes > 4`:

| id_examples (10 of 12)                | `nodes` | `node4` | expected |
|---------------------------------------|---------|---------|----------|
| see `audit_summary.csv` row for count | 1, 3, 4 | 1       | 0 (false positive) |
| see `audit_summary.csv` row for count | 5, 8, 9 | 0       | 1 (false negative) |

These are residual data-entry errors from the original NCCTG case-report
forms. They were *not* corrected when the dataset was deposited in the R
`survival` package.

**How we handle it.**

- Primary analyses use `nodes` (the continuous count).
- Where a binary node burden is needed, we recompute `nodes_high =
  (nodes > 4).astype(int)` rather than using the column-level `node4`.
- The Q1 subgroup forest uses `nodes_high`.
- The 12 contradictory rows remain in the analytic frame and are reported
  alongside the audit summary.

### W2. Follow-up window

The original Moertel 1990 paper reported ~3-year follow-up. The dataset in
`survival`/`ForCausality` reflects the extended 5–9 year follow-up from the
1995 *Annals* update. We use the full follow-up for Q1 (Cox + RMST(5y)),
which gives the *same* HR ≈ 0.67 to two decimal places as the 1990 paper
under PH.

This is *important for Methods*: do not claim "we reproduce the 1990
result" if the follow-up window is different. The honest claim is "we
reproduce the 1990 HR using the 1995 follow-up extension, which preserves
the result under PH."

### W3. Missingness

`nodes` (18 missing, 1.9%) and `differ` (23 missing, 2.5%) are the only
covariates with any missing values. Both come from the pathology report.
Likely missing-completely-at-random with respect to treatment assignment,
because pathologists could not know `rx` when generating the report.

**Primary analysis** uses per-arm median imputation. **Sensitivity**
(reported in §S2) uses complete-case analysis.

### W4. `study` column

Always `1`. Harmless legacy field. We drop it from analytic frames.

---

## Reproducibility

To regenerate this CSV from scratch (requires R with the `ForCausality`
package installed):

```bash
python scripts/pull_colon_data.py
```

The script verifies row count, unique-ID count, etype split, and `rx`
levels before writing. If any verification gate fails, the script exits
non-zero and does *not* overwrite `data/colon.csv`.

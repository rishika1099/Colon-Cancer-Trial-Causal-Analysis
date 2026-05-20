# IEEE Conference Report

This directory contains the IEEE conference-format write-up of the causal
re-analysis.

| File | Description |
|---|---|
| `ieee_report.tex` | Source — IEEEtran two-column conference template |
| `ieee_report.pdf` | Compiled PDF (5 pages, US Letter) |

## Build

Either toolchain works:

```bash
# Tectonic (recommended — self-contained, fetches packages on demand)
tectonic ieee_report.tex

# Or TeX Live
pdflatex ieee_report.tex
pdflatex ieee_report.tex   # second pass for cross-refs
```

The references are inlined in the `.tex` (no separate `.bib` step) for
build simplicity. A BibTeX version is available in
`../manuscript/references.bib` if you prefer to convert.

## Contents at a glance

| Section | What it covers |
|---|---|
| Abstract + Keywords | Five-question framing; headline numbers |
| I. Introduction | Problem statement, contributions |
| II. Data and Preliminaries | Audit, notation, two DAGs |
| III. Q1 — Randomized ATE | KM, Cox, RMST, PH check, anchor pass |
| IV. Q2 — Back-door + bad-control | Five estimators converge; bad-control flips |
| V. Q3 — Heterogeneous effects | Meta-learners + causal forest + CATE-by-nodes |
| VI. Q4 — Mediation | Imai-Keele-Tingley + ρ sensitivity |
| VII. Q5 — Transportability | Cole-Stuart IOSW + Dahabreh bound |
| VIII. Sensitivity Synthesis | VanderWeele-Ding E-values |
| IX. Discussion | Three pedagogical findings |
| X. Conclusion + Acknowledgment + References | — |

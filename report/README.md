# IEEE Conference Report

This directory contains the IEEE conference-format write-up of the
project.

| File | Description |
|---|---|
| `ieee_report.tex` | Source (IEEEtran two-column conference template). |
| `ieee_report.pdf` | Compiled PDF, 5 pages, US Letter. |

## Build

Either toolchain works:

```bash
# Tectonic (recommended: self-contained, no system TeX install required)
tectonic ieee_report.tex

# Or TeX Live
pdflatex ieee_report.tex
pdflatex ieee_report.tex   # second pass for cross-refs
```

The references are inlined in the `.tex` source, so no separate BibTeX
step is needed. A BibTeX version of the same references is available in
`../manuscript/references.bib`.

## Contents

| Section | Topic |
|---|---|
| Abstract + keywords | Five-question framing and headline numbers |
| I. Introduction | Background, contributions |
| II. Data and Preliminaries | Audit, notation, DAGs |
| III. Q1 — Randomized ATE | KM, Cox, RMST, PH check |
| IV. Q2 — Back-door + bad-control | Five estimators plus the bad-control comparison |
| V. Q3 — Heterogeneous effects | Meta-learners, causal forest, CATE-by-nodes |
| VI. Q4 — Mediation | Imai-Keele-Tingley + ρ-sensitivity |
| VII. Q5 — Transportability | Cole-Stuart IOSW + Dahabreh bound |
| VIII. Sensitivity Synthesis | VanderWeele-Ding E-values |
| IX. Discussion | Three substantive findings |
| X. Conclusion + References | — |

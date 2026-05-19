"""Pull ForCausality::Colon_df via rpy2 and write data/colon.csv.

Source of truth: the R package `ForCausality` (Toby Codigos), which is a
curated copy of `survival::colon` — the Moertel et al. NEJM 1990 trial.

Verification gates (all must pass):
    * 1,858 rows (2 per patient: etype=1 recurrence, etype=2 death).
    *   929 unique patient IDs.
    *   929 rows with etype == 2 (one death-row record per patient).
    * `rx` levels = {Obs, Lev, Lev+5FU}.

Run from project root:
    python scripts/pull_colon_data.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "colon.csv"

EXPECTED_ROWS = 1858
EXPECTED_IDS = 929
EXPECTED_DEATHS = 929
EXPECTED_RX = {"Obs", "Lev", "Lev+5FU"}


def pull_via_rpy2() -> "pandas.DataFrame":
    """Pull Colon_df through rpy2."""
    from rpy2.robjects import r, pandas2ri
    from rpy2.robjects.packages import importr

    pandas2ri.activate()
    importr("ForCausality")
    # The dataset is lazy-loaded; force materialization.
    df = pandas2ri.rpy2py(r("ForCausality::Colon_df"))
    return df


def pull_via_subprocess() -> "pandas.DataFrame":
    """Fallback: shell out to Rscript and write a temp CSV.

    Used when rpy2 is not installed yet (e.g. fresh clone before
    `pip install -r requirements.txt`).
    """
    import subprocess
    import tempfile
    import pandas as pd

    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
        tmp_path = tmp.name
    rscript = (
        "suppressPackageStartupMessages(library(ForCausality));"
        f"write.csv(Colon_df, '{tmp_path}', row.names = FALSE)"
    )
    subprocess.run(["Rscript", "-e", rscript], check=True)
    df = pd.read_csv(tmp_path)
    Path(tmp_path).unlink(missing_ok=True)
    return df


def verify(df) -> None:
    rows = len(df)
    ids = df["id"].nunique()
    deaths = int((df["etype"] == 2).sum())
    rx = set(df["rx"].astype(str).unique())

    print(f"  rows                : {rows}    (expected {EXPECTED_ROWS})")
    print(f"  unique patient IDs  : {ids}     (expected {EXPECTED_IDS})")
    print(f"  etype==2 records    : {deaths}  (expected {EXPECTED_DEATHS})")
    print(f"  rx levels           : {sorted(rx)}")

    if rows != EXPECTED_ROWS:
        sys.exit(f"FAIL: expected {EXPECTED_ROWS} rows, got {rows}")
    if ids != EXPECTED_IDS:
        sys.exit(f"FAIL: expected {EXPECTED_IDS} unique IDs, got {ids}")
    if deaths != EXPECTED_DEATHS:
        sys.exit(f"FAIL: expected {EXPECTED_DEATHS} etype==2 rows, got {deaths}")
    if rx != EXPECTED_RX:
        sys.exit(f"FAIL: expected rx={sorted(EXPECTED_RX)}, got {sorted(rx)}")
    # Spot-check expected schema.
    expected_cols = {
        "id", "study", "rx", "sex", "age", "obstruct", "perfor", "adhere",
        "nodes", "status", "differ", "extent", "surg", "node4", "time", "etype",
    }
    missing = expected_cols - set(df.columns)
    if missing:
        sys.exit(f"FAIL: missing columns: {sorted(missing)}")
    print("  schema              : ok")


def main() -> None:
    print("Pulling ForCausality::Colon_df ...")
    try:
        df = pull_via_rpy2()
        print("  source : rpy2")
    except Exception as exc:  # noqa: BLE001
        print(f"  rpy2 unavailable ({exc.__class__.__name__}); falling back to Rscript.")
        df = pull_via_subprocess()
        print("  source : Rscript subprocess")

    # rpy2 sometimes returns rx as a categorical with integer codes — coerce
    # to the human-readable factor labels.
    if df["rx"].dtype.kind in {"i", "u"}:
        df["rx"] = (
            df["rx"].map({1: "Obs", 2: "Lev", 3: "Lev+5FU"})
        )
    df["rx"] = df["rx"].astype(str)

    print("Verifying:")
    verify(df)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    print(f"\nWrote {OUT.relative_to(ROOT)}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()

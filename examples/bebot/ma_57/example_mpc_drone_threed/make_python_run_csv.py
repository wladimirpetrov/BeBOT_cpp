#!/usr/bin/env python3
"""
Build python_run.csv from the current SOS/MATLAB wall solution plus an existing
mpc_vectors_from_smoothed_trajectory.csv file.

The SOS part is read from wall_sos_la3_generated_solution.csv produced by
run_sos_codegen_wall_ocp_3d.m. For the SOS rows:
  - input_x/input_y/input_z come from ux1/uy1/uz1
  - init_px/init_py/init_pz come from x1/y1/z1
  - init_psi is set to 0
  - all px_c*, py_c*, pz_c*, psi_c* coefficient columns are set to 0

Then the existing MPC-vector rows are appended after the SOS rows, with point_id
renumbered continuously and arc_length offset so it continues from the SOS path.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


BASE_COLUMNS = [
    "point_id",
    "arc_length",
    "input_x",
    "input_y",
    "input_z",
    "init_px",
    "init_py",
    "init_pz",
    "init_psi",
]

COEFF_COLUMNS = (
    [f"px_c{i}" for i in range(5)]
    + [f"py_c{i}" for i in range(5)]
    + [f"pz_c{i}" for i in range(5)]
    + [f"psi_c{i}" for i in range(5)]
)

OUTPUT_COLUMNS = BASE_COLUMNS + COEFF_COLUMNS


def parse_args() -> argparse.Namespace:
    default_sos = Path.home() / "dev" / "SOS" / "fw_sos_software" / "my_problems" / "wall_sos_la3_codegen" / "wall_sos_la3_generated_solution.csv"

    parser = argparse.ArgumentParser(
        description="Create python_run.csv from current SOS solution + existing MPC-vector CSV."
    )
    parser.add_argument(
        "--sos-csv",
        default=str(default_sos),
        help="Path to SOS/MATLAB generated solution CSV. Default: %(default)s",
    )
    parser.add_argument(
        "--mpc-csv",
        default="mpc_vectors_from_smoothed_trajectory.csv",
        help="Path to existing mpc_vectors_from_smoothed_trajectory.csv. Default: %(default)s",
    )
    parser.add_argument(
        "--output",
        default="python_run.csv",
        help="Output CSV path. Default: %(default)s",
    )
    return parser.parse_args()


def require_columns(df: pd.DataFrame, cols: Iterable[str], source_name: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"{source_name} is missing required columns {missing}.\n"
            f"Available columns are: {list(df.columns)}"
        )


def cumulative_arc_length(x: pd.Series, y: pd.Series, z: pd.Series) -> np.ndarray:
    x_arr = x.to_numpy(dtype=float)
    y_arr = y.to_numpy(dtype=float)
    z_arr = z.to_numpy(dtype=float)

    if len(x_arr) == 0:
        return np.array([], dtype=float)

    ds = np.sqrt(np.diff(x_arr) ** 2 + np.diff(y_arr) ** 2 + np.diff(z_arr) ** 2)
    return np.concatenate(([0.0], np.cumsum(ds)))


def enforce_output_columns(df: pd.DataFrame, fill_missing_with_zero: bool = True) -> pd.DataFrame:
    out = df.copy()

    for col in OUTPUT_COLUMNS:
        if col not in out.columns:
            if fill_missing_with_zero:
                out[col] = 0.0
            else:
                raise ValueError(f"Missing required output column: {col}")

    return out[OUTPUT_COLUMNS]


def build_sos_output_rows(sos_csv: Path) -> pd.DataFrame:
    sos_df = pd.read_csv(sos_csv)

    require_columns(
        sos_df,
        ["x1", "y1", "z1", "ux1", "uy1", "uz1"],
        str(sos_csv),
    )

    n = len(sos_df)
    out = pd.DataFrame()

    # Use 0-based point ids. Change to np.arange(1, n + 1) if your ROS player expects 1-based ids.
    out["point_id"] = np.arange(n, dtype=int)

    out["arc_length"] = cumulative_arc_length(sos_df["x1"], sos_df["y1"], sos_df["z1"])

    out["input_x"] = sos_df["ux1"].astype(float)
    out["input_y"] = sos_df["uy1"].astype(float)
    out["input_z"] = sos_df["uz1"].astype(float)

    out["init_px"] = sos_df["x1"].astype(float)
    out["init_py"] = sos_df["y1"].astype(float)
    out["init_pz"] = sos_df["z1"].astype(float)
    out["init_psi"] = 0.0

    for col in COEFF_COLUMNS:
        out[col] = 0.0

    return enforce_output_columns(out)


def append_mpc_rows_after_sos(sos_out: pd.DataFrame, mpc_csv: Path) -> pd.DataFrame:
    mpc_df = pd.read_csv(mpc_csv)

    # Keep existing MPC vector columns where present. Missing required columns are filled with zero.
    mpc_out = enforce_output_columns(mpc_df, fill_missing_with_zero=True)

    n_sos = len(sos_out)
    n_mpc = len(mpc_out)

    # Continue point_id after SOS rows.
    mpc_out["point_id"] = np.arange(n_sos, n_sos + n_mpc, dtype=int)

    # Continue arc_length after the SOS arc length while preserving relative MPC arc lengths.
    if n_sos > 0 and n_mpc > 0:
        sos_last_arc = float(sos_out["arc_length"].iloc[-1])
        mpc_first_arc = float(mpc_out["arc_length"].iloc[0])
        mpc_out["arc_length"] = sos_last_arc + (mpc_out["arc_length"].astype(float) - mpc_first_arc)

    combined = pd.concat([sos_out, mpc_out], ignore_index=True)
    return enforce_output_columns(combined)


def main() -> None:
    args = parse_args()

    sos_csv = Path(args.sos_csv).expanduser().resolve()
    mpc_csv = Path(args.mpc_csv).expanduser().resolve()
    output_csv = Path(args.output).expanduser().resolve()

    if not sos_csv.is_file():
        raise FileNotFoundError(f"SOS solution CSV not found: {sos_csv}")
    if not mpc_csv.is_file():
        raise FileNotFoundError(f"MPC vectors CSV not found: {mpc_csv}")

    sos_out = build_sos_output_rows(sos_csv)
    combined = append_mpc_rows_after_sos(sos_out, mpc_csv)

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_csv, index=False)

    print(f"Saved: {output_csv}")
    print(f"SOS rows:      {len(sos_out)}")
    print(f"MPC rows:      {len(combined) - len(sos_out)}")
    print(f"Total rows:    {len(combined)}")
    print(f"First point_id: {combined['point_id'].iloc[0] if len(combined) else 'N/A'}")
    print(f"Last point_id:  {combined['point_id'].iloc[-1] if len(combined) else 'N/A'}")


if __name__ == "__main__":
    main()

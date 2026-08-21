#!/usr/bin/env python3

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DATASET = "single_drone_single_obs_training_dataset.csv"
CURVE_SAMPLES = 250


def next_available_path(filename):
    path = Path(filename)

    if not path.exists():
        return path

    stem = path.stem
    suffix = path.suffix

    for i in range(1, 10000):
        candidate = path.with_name(f"{stem}_{i:03d}{suffix}")
        if not candidate.exists():
            return candidate

    raise RuntimeError(f"Could not find free output filename for {filename}")


def bernstein_curve(control_points, number_of_samples=CURVE_SAMPLES):
    """
    Geometric Bernstein curve evaluation on s in [0,1].

    This is enough for checking the geometric path produced by the saved
    Bernstein coefficients. It does not assume anything about physical time.
    """
    cp = np.asarray(control_points, dtype=float)
    degree = len(cp) - 1

    s = np.linspace(0.0, 1.0, number_of_samples)

    curve = np.zeros_like(s)

    for i in range(degree + 1):
        basis = (
            math.comb(degree, i)
            * (s ** i)
            * ((1.0 - s) ** (degree - i))
        )

        curve += cp[i] * basis

    return s, curve


def get_control_points(row, prefix, degree):
    return np.array(
        [row[f"{prefix}_cp_{i}"] for i in range(degree + 1)],
        dtype=float,
    )


def trajectory_from_row(row):
    degree = int(row["N"])

    _, x = bernstein_curve(
        get_control_points(row, "px", degree)
    )

    _, y = bernstein_curve(
        get_control_points(row, "py", degree)
    )

    _, z = bernstein_curve(
        get_control_points(row, "pz", degree)
    )

    _, psi = bernstein_curve(
        get_control_points(row, "psi", degree)
    )

    return x, y, z, psi


def minimum_clearance(row):
    x, y, z, _ = trajectory_from_row(row)

    dx = x - row["sphere_x"]
    dy = y - row["sphere_y"]
    dz = z - row["sphere_z"]

    distance = np.sqrt(
        dx * dx
        + dy * dy
        + dz * dz
    )

    min_distance = float(np.min(distance))
    clearance = min_distance - float(row["sphere_radius"])

    return min_distance, clearance


def print_dataset_summary(df):
    print("============================================================")
    print("DATASET SUMMARY")
    print("============================================================")
    print(f"Rows / successful OCP solutions : {len(df)}")

    if "sample_id" in df.columns:
        print(
            f"Unique sample_id values         : "
            f"{df['sample_id'].nunique()}"
        )

    print("\nCase counts:")
    print(df["case_type"].value_counts().to_string())

    print("\nCase percentages:")
    percentages = (
        100.0
        * df["case_type"].value_counts(normalize=True)
    )
    print(percentages.round(2).to_string())

    print("\nObstacle radius:")
    print(
        f"  min  = {df['sphere_radius'].min():.6f}\n"
        f"  mean = {df['sphere_radius'].mean():.6f}\n"
        f"  max  = {df['sphere_radius'].max():.6f}"
    )

    rho_reconstructed = np.sqrt(
        df["sphere_y"].to_numpy() ** 2
        + df["sphere_z"].to_numpy() ** 2
    )

    rho_error = np.max(
        np.abs(
            rho_reconstructed
            - df["sphere_rho"].to_numpy()
        )
    )

    print(
        "\nMaximum |saved rho - sqrt(y^2+z^2)|: "
        f"{rho_error:.3e}"
    )

    print("============================================================\n")


def plot_obstacle_distribution(df):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    for case_type in ["blocking", "near_miss", "clear"]:
        subset = df[df["case_type"] == case_type]

        ax.scatter(
            subset["sphere_x"],
            subset["sphere_y"],
            subset["sphere_z"],
            s=42,
            alpha=0.65,
            label=f"{case_type} ({len(subset)})",
        )

    ax.plot(
        [0.0, 1.0],
        [0.0, 0.0],
        [0.0, 0.0],
        linewidth=2.5,
        label="Direct start-target line",
    )

    ax.scatter(
        [0.0],
        [0.0],
        [0.0],
        marker="*",
        s=150,
        label="Start",
    )

    ax.scatter(
        [1.0],
        [0.0],
        [0.0],
        marker="X",
        s=120,
        label="Target",
    )

    ax.set_xlabel("Canonical x")
    ax.set_ylabel("Canonical y")
    ax.set_zlabel("Canonical z")

    ax.set_title(
        "Actual saved obstacle distribution"
    )

    ax.legend(fontsize=8)

    fig.tight_layout()

    output = next_available_path(
        "single_obs_dataset_obstacles_3d.png"
    )

    fig.savefig(
        output,
        dpi=220,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Saved: {output}")


def plot_sampling_regions(df):
    fig, ax = plt.subplots(figsize=(10, 7))

    for case_type in ["blocking", "near_miss", "clear"]:
        subset = df[df["case_type"] == case_type]

        ratio = (
            subset["sphere_rho"]
            / subset["sphere_radius"]
        )

        ax.scatter(
            subset["sphere_x"],
            ratio,
            s=42,
            alpha=0.65,
            label=case_type,
        )

    ax.axhline(
        1.0,
        linestyle="--",
        linewidth=2.0,
        label="rho/r = 1 collision boundary",
    )

    ax.set_xlabel(
        "Obstacle x-position along start-target segment"
    )

    ax.set_ylabel(
        "rho / sphere radius"
    )

    ax.set_title(
        "Actual saved geometric case distribution"
    )

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 3.2)

    ax.grid(alpha=0.25)
    ax.legend()

    fig.tight_layout()

    output = next_available_path(
        "single_obs_dataset_case_regions.png"
    )

    fig.savefig(
        output,
        dpi=220,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Saved: {output}")


def choose_representative_case(df, case_type):
    subset = df[df["case_type"] == case_type].copy()

    if len(subset) == 0:
        return None

    ratio = (
        subset["sphere_rho"]
        / subset["sphere_radius"]
    )

    ratio_mid = float(ratio.median())

    # Prefer an obstacle near the middle of the start-target segment
    # and near the median geometry of its category.
    score = (
        np.abs(subset["sphere_x"] - 0.5)
        + 0.25 * np.abs(ratio - ratio_mid)
    )

    index = score.idxmin()

    return subset.loc[index]


def plot_representative_case(row):
    case_type = str(row["case_type"])

    x, y, z, _ = trajectory_from_row(row)

    sphere_x = float(row["sphere_x"])
    sphere_y = float(row["sphere_y"])
    sphere_z = float(row["sphere_z"])
    radius = float(row["sphere_radius"])

    min_distance, clearance = minimum_clearance(row)

    u = np.linspace(
        0.0,
        2.0 * np.pi,
        45,
    )

    v = np.linspace(
        0.0,
        np.pi,
        28,
    )

    sphere_surface_x = (
        sphere_x
        + radius
        * np.outer(
            np.cos(u),
            np.sin(v),
        )
    )

    sphere_surface_y = (
        sphere_y
        + radius
        * np.outer(
            np.sin(u),
            np.sin(v),
        )
    )

    sphere_surface_z = (
        sphere_z
        + radius
        * np.outer(
            np.ones_like(u),
            np.cos(v),
        )
    )

    fig = plt.figure(figsize=(10, 8))

    ax = fig.add_subplot(
        111,
        projection="3d",
    )

    ax.plot(
        x,
        y,
        z,
        linewidth=2.8,
        label="Optimized Bernstein path",
    )

    ax.plot(
        [0.0, 1.0],
        [0.0, 0.0],
        [0.0, 0.0],
        linestyle="--",
        linewidth=1.7,
        label="Straight start-target path",
    )

    ax.plot_wireframe(
        sphere_surface_x,
        sphere_surface_y,
        sphere_surface_z,
        rstride=2,
        cstride=2,
        linewidth=0.55,
        alpha=0.6,
    )

    ax.scatter(
        [0.0],
        [0.0],
        [0.0],
        marker="*",
        s=150,
        label="Start",
    )

    ax.scatter(
        [1.0],
        [0.0],
        [0.0],
        marker="X",
        s=120,
        label="Target",
    )

    ax.scatter(
        [sphere_x],
        [sphere_y],
        [sphere_z],
        s=75,
        label="Sphere center",
    )

    ax.set_xlabel("Canonical x")
    ax.set_ylabel("Canonical y")
    ax.set_zlabel("Canonical z")

    sample_id = int(row["sample_id"])

    ax.set_title(
        f"{case_type}: sample {sample_id} | "
        f"clearance = {clearance:.4f}"
    )

    ax.legend(fontsize=8)

    fig.tight_layout()

    output = next_available_path(
        f"single_obs_dataset_example_{case_type}.png"
    )

    fig.savefig(
        output,
        dpi=220,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(
        f"Saved: {output}\n"
        f"  sample_id       = {sample_id}\n"
        f"  sphere center   = "
        f"[{sphere_x:.4f}, {sphere_y:.4f}, {sphere_z:.4f}]\n"
        f"  sphere radius   = {radius:.4f}\n"
        f"  min distance    = {min_distance:.4f}\n"
        f"  dense clearance = {clearance:.4f}\n"
    )


def plot_clearance_distribution(df):
    records = []

    print(
        "Computing dense geometric clearance "
        "for every saved trajectory..."
    )

    for _, row in df.iterrows():
        _, clearance = minimum_clearance(row)

        records.append(
            {
                "case_type": row["case_type"],
                "sphere_x": row["sphere_x"],
                "clearance": clearance,
            }
        )

    clearance_df = pd.DataFrame(records)

    fig, ax = plt.subplots(figsize=(10, 7))

    for case_type in ["blocking", "near_miss", "clear"]:
        subset = clearance_df[
            clearance_df["case_type"]
            == case_type
        ]

        ax.scatter(
            subset["sphere_x"],
            subset["clearance"],
            s=42,
            alpha=0.65,
            label=case_type,
        )

    ax.axhline(
        0.0,
        linestyle="--",
        linewidth=2.0,
        label="Collision boundary",
    )

    ax.set_xlabel(
        "Obstacle x-position"
    )

    ax.set_ylabel(
        "Dense path clearance to sphere"
    )

    ax.set_title(
        "Collision check for all saved optimized trajectories"
    )

    ax.grid(alpha=0.25)
    ax.legend()

    fig.tight_layout()

    output = next_available_path(
        "single_obs_dataset_clearance.png"
    )

    fig.savefig(
        output,
        dpi=220,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Saved: {output}")

    print("\nDense trajectory-clearance summary:")
    print(
        clearance_df
        .groupby("case_type")["clearance"]
        .agg(["min", "mean", "max"])
        .to_string()
    )

    collisions = clearance_df[
        clearance_df["clearance"]
        < -5e-3
    ]

    print(
        "\nTrajectories with penetration > 0.005: "
        f"{len(collisions)}"
    )

    return clearance_df


def main():
    dataset_path = Path(DATASET)

    if not dataset_path.exists():
        raise FileNotFoundError(
            f"Could not find {DATASET} in the current directory."
        )

    df = pd.read_csv(dataset_path)

    required_columns = [
        "sample_id",
        "case_type",
        "N",
        "start_x",
        "start_y",
        "start_z",
        "target_x",
        "target_y",
        "target_z",
        "sphere_x",
        "sphere_y",
        "sphere_z",
        "sphere_radius",
        "sphere_rho",
    ]

    for column in required_columns:
        if column not in df.columns:
            raise ValueError(
                f"Dataset is missing required column: {column}"
            )

    print_dataset_summary(df)

    plot_obstacle_distribution(df)

    plot_sampling_regions(df)

    clearance_df = plot_clearance_distribution(df)

    for case_type in [
        "blocking",
        "near_miss",
        "clear",
    ]:
        row = choose_representative_case(
            df,
            case_type,
        )

        if row is not None:
            plot_representative_case(row)

    print(
        "\nDone. Inspect the generated PNG files "
        "in the current directory."
    )


if __name__ == "__main__":
    main()

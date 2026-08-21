#!/usr/bin/env python3

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ============================================================
# INPUT
# Produced by compare_nlp_vs_nn.py
# ============================================================
CONTROL_POINT_CSV = "nlp_vs_nn_control_points.csv"

# Same test case used by compare_nlp_vs_nn.py.
# Change these four numbers if you run the comparison with
# different --sphere-x, --sphere-y, --sphere-z, --radius values.
SPHERE_X = 0.55
SPHERE_Y = 0.08
SPHERE_Z = -0.05
SPHERE_RADIUS = 0.14

START = np.array([0.0, 0.0, 0.0])
TARGET = np.array([1.0, 0.0, 0.0])

CURVE_SAMPLES = 500


# ============================================================
# Bernstein evaluation
# ============================================================
def bernstein_curve(control_points, samples=CURVE_SAMPLES):
    control_points = np.asarray(
        control_points,
        dtype=float
    )

    N = len(control_points) - 1

    s = np.linspace(
        0.0,
        1.0,
        samples
    )

    curve = np.zeros_like(s)

    for i in range(N + 1):
        basis = (
            math.comb(N, i)
            * (s ** i)
            * ((1.0 - s) ** (N - i))
        )

        curve += (
            control_points[i]
            * basis
        )

    return curve


# ============================================================
# Load one variable from the comparison CSV
# ============================================================
def get_control_points(df, variable, column):
    subset = (
        df[df["variable"] == variable]
        .sort_values("cp_index")
    )

    if len(subset) == 0:
        raise ValueError(
            f"No rows found for variable '{variable}'."
        )

    return subset[column].to_numpy(
        dtype=float
    )


# ============================================================
# Dense clearance
# ============================================================
def compute_clearance(x, y, z):
    distance = np.sqrt(
        (x - SPHERE_X) ** 2
        + (y - SPHERE_Y) ** 2
        + (z - SPHERE_Z) ** 2
    )

    minimum_distance = float(
        np.min(distance)
    )

    clearance = (
        minimum_distance
        - SPHERE_RADIUS
    )

    return minimum_distance, clearance


# ============================================================
# Main
# ============================================================
def main():
    csv_path = Path(
        CONTROL_POINT_CSV
    )

    if not csv_path.exists():
        raise FileNotFoundError(
            f"Could not find {CONTROL_POINT_CSV}"
        )

    df = pd.read_csv(
        csv_path
    )

    required_columns = [
        "variable",
        "cp_index",
        "nlp",
        "nn",
    ]

    for column in required_columns:
        if column not in df.columns:
            raise ValueError(
                f"Missing required column: {column}"
            )

    # --------------------------------------------------------
    # NLP Bernstein control points
    # --------------------------------------------------------
    nlp_px = get_control_points(
        df,
        "px",
        "nlp"
    )

    nlp_py = get_control_points(
        df,
        "py",
        "nlp"
    )

    nlp_pz = get_control_points(
        df,
        "pz",
        "nlp"
    )

    # --------------------------------------------------------
    # NN Bernstein control points
    # --------------------------------------------------------
    nn_px = get_control_points(
        df,
        "px",
        "nn"
    )

    nn_py = get_control_points(
        df,
        "py",
        "nn"
    )

    nn_pz = get_control_points(
        df,
        "pz",
        "nn"
    )

    # --------------------------------------------------------
    # Reconstruct trajectories
    # --------------------------------------------------------
    nlp_x = bernstein_curve(nlp_px)
    nlp_y = bernstein_curve(nlp_py)
    nlp_z = bernstein_curve(nlp_pz)

    nn_x = bernstein_curve(nn_px)
    nn_y = bernstein_curve(nn_py)
    nn_z = bernstein_curve(nn_pz)

    # --------------------------------------------------------
    # Clearance
    # --------------------------------------------------------
    nlp_min_distance, nlp_clearance = (
        compute_clearance(
            nlp_x,
            nlp_y,
            nlp_z
        )
    )

    nn_min_distance, nn_clearance = (
        compute_clearance(
            nn_x,
            nn_y,
            nn_z
        )
    )

    print(
        "============================================="
    )

    print(
        "NLP vs NN trajectory comparison"
    )

    print(
        "============================================="
    )

    print(
        f"NLP minimum sphere distance : "
        f"{nlp_min_distance:.10f}"
    )

    print(
        f"NLP clearance               : "
        f"{nlp_clearance:.10f}"
    )

    print(
        f"NN minimum sphere distance  : "
        f"{nn_min_distance:.10f}"
    )

    print(
        f"NN clearance                : "
        f"{nn_clearance:.10f}"
    )

    # --------------------------------------------------------
    # Sphere surface
    # --------------------------------------------------------
    u = np.linspace(
        0.0,
        2.0 * np.pi,
        50
    )

    v = np.linspace(
        0.0,
        np.pi,
        30
    )

    sphere_surface_x = (
        SPHERE_X
        + SPHERE_RADIUS
        * np.outer(
            np.cos(u),
            np.sin(v)
        )
    )

    sphere_surface_y = (
        SPHERE_Y
        + SPHERE_RADIUS
        * np.outer(
            np.sin(u),
            np.sin(v)
        )
    )

    sphere_surface_z = (
        SPHERE_Z
        + SPHERE_RADIUS
        * np.outer(
            np.ones_like(u),
            np.cos(v)
        )
    )

    # ========================================================
    # Plot BOTH trajectories on the SAME figure
    # ========================================================
    fig = plt.figure(
        figsize=(10, 8)
    )

    ax = fig.add_subplot(
        111,
        projection="3d"
    )

    # NLP
    ax.plot(
        nlp_x,
        nlp_y,
        nlp_z,
        linewidth=3.0,
        label=(
            "C++ NLP "
            f"(clearance={nlp_clearance:.4f})"
        )
    )

    # NN
    ax.plot(
        nn_x,
        nn_y,
        nn_z,
        linewidth=3.0,
        linestyle="--",
        label=(
            "Neural network "
            f"(clearance={nn_clearance:.4f})"
        )
    )

    # Direct start-target line
    ax.plot(
        [START[0], TARGET[0]],
        [START[1], TARGET[1]],
        [START[2], TARGET[2]],
        linestyle=":",
        linewidth=1.5,
        label="Direct start-target line"
    )

    # Sphere
    ax.plot_wireframe(
        sphere_surface_x,
        sphere_surface_y,
        sphere_surface_z,
        rstride=2,
        cstride=2,
        linewidth=0.55,
        alpha=0.55
    )

    # Start
    ax.scatter(
        [START[0]],
        [START[1]],
        [START[2]],
        s=150,
        marker="*",
        label="Start"
    )

    # Target
    ax.scatter(
        [TARGET[0]],
        [TARGET[1]],
        [TARGET[2]],
        s=120,
        marker="X",
        label="Target"
    )

    # Sphere center
    ax.scatter(
        [SPHERE_X],
        [SPHERE_Y],
        [SPHERE_Z],
        s=75,
        marker="o",
        label="Sphere center"
    )

    ax.set_xlabel(
        "Canonical x"
    )

    ax.set_ylabel(
        "Canonical y"
    )

    ax.set_zlabel(
        "Canonical z"
    )

    ax.set_title(
        "C++ NLP vs Neural Network\n"
        "Bernstein trajectory comparison"
    )

    ax.legend(
        fontsize=8
    )

    # Same useful limits for this canonical problem.
    ax.set_xlim(
        -0.05,
        1.05
    )

    # Automatically choose y/z limits with some margin so that
    # both trajectories and the obstacle are visible.
    all_y = np.concatenate(
        [
            nlp_y,
            nn_y,
            np.array(
                [
                    SPHERE_Y - SPHERE_RADIUS,
                    SPHERE_Y + SPHERE_RADIUS
                ]
            )
        ]
    )

    all_z = np.concatenate(
        [
            nlp_z,
            nn_z,
            np.array(
                [
                    SPHERE_Z - SPHERE_RADIUS,
                    SPHERE_Z + SPHERE_RADIUS
                ]
            )
        ]
    )

    y_min = float(all_y.min())
    y_max = float(all_y.max())

    z_min = float(all_z.min())
    z_max = float(all_z.max())

    y_margin = max(
        0.05,
        0.15 * (y_max - y_min)
    )

    z_margin = max(
        0.05,
        0.15 * (z_max - z_min)
    )

    ax.set_ylim(
        y_min - y_margin,
        y_max + y_margin
    )

    ax.set_zlim(
        z_min - z_margin,
        z_max + z_margin
    )

    fig.tight_layout()

    output_png = (
        "nlp_vs_nn_trajectory_plot.png"
    )

    fig.savefig(
        output_png,
        dpi=220,
        bbox_inches="tight"
    )

    print(
        f"\nSaved plot: {output_png}"
    )

    # Show interactively.
    plt.show()


if __name__ == "__main__":
    main()
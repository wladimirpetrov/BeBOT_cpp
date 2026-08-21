#!/usr/bin/env python3

import argparse
import ctypes
import math
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from torch import nn


# ============================================================
# Neural-network architecture
# Must match train_single_obs_mlp.py exactly.
# ============================================================
class ResidualBlock(nn.Module):
    def __init__(self, width):
        super().__init__()

        self.net = nn.Sequential(
            nn.Linear(width, width),
            nn.SiLU(),
            nn.Linear(width, width),
        )

        self.activation = nn.SiLU()

    def forward(self, x):
        return self.activation(
            x + self.net(x)
        )


class BernsteinMLP(nn.Module):
    def __init__(
        self,
        input_dim,
        output_dim,
        width=128,
        residual_blocks=3,
    ):
        super().__init__()

        layers = [
            nn.Linear(input_dim, 64),
            nn.SiLU(),
            nn.Linear(64, width),
            nn.SiLU(),
        ]

        for _ in range(residual_blocks):
            layers.append(
                ResidualBlock(width)
            )

        layers.extend(
            [
                nn.Linear(width, 64),
                nn.SiLU(),
                nn.Linear(64, output_dim),
            ]
        )

        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)


# ============================================================
# File helper: do not overwrite previous comparison results.
# ============================================================
def non_overwriting_path(filename):
    path = Path(filename)

    if not path.exists():
        return path

    for i in range(1, 10000):
        candidate = path.with_name(
            "{}_{:03d}{}".format(
                path.stem,
                i,
                path.suffix,
            )
        )

        if not candidate.exists():
            return candidate

    raise RuntimeError(
        "Could not find free filename for {}".format(
            filename
        )
    )


# ============================================================
# Load C++ NLP shared library.
# ============================================================
def load_nlp_library(library_path):
    library_path = Path(library_path).resolve()

    if not library_path.exists():
        raise FileNotFoundError(
            "Could not find NLP shared library: {}".format(
                library_path
            )
        )

    lib = ctypes.CDLL(str(library_path))

    # create_point_set_problem:
    #   1 int + 29 doubles
    lib.create_point_set_problem.argtypes = (
        [ctypes.c_int]
        + [ctypes.c_double] * 29
    )

    lib.create_point_set_problem.restype = (
        ctypes.c_void_p
    )

    lib.solve_point_set_problem.argtypes = [
        ctypes.c_void_p
    ]

    lib.solve_point_set_problem.restype = (
        ctypes.c_int
    )

    lib.get_solution_size.argtypes = [
        ctypes.c_void_p
    ]

    lib.get_solution_size.restype = (
        ctypes.c_int
    )

    lib.get_solution.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
    ]

    lib.get_solution.restype = (
        ctypes.c_int
    )

    lib.get_final_objective_value.argtypes = [
        ctypes.c_void_p
    ]

    lib.get_final_objective_value.restype = (
        ctypes.c_double
    )

    lib.destroy_point_set_problem.argtypes = [
        ctypes.c_void_p
    ]

    lib.destroy_point_set_problem.restype = None

    return lib


# ============================================================
# Solve canonical one-sphere OCP using the C++ NLP.
#
# Canonical problem:
#   start  = [0,0,0]
#   target = [1,0,0]
#   initial velocity = 0
#   initial yaw = target yaw = 0
# ============================================================
def solve_cpp_nlp(
    lib,
    sphere_x,
    sphere_y,
    sphere_z,
    sphere_radius,
    N=4,
    tf=0.9,
):
    inf = float("inf")

    # Canonical bounds used by the data-collection caller.
    px_min = -inf
    px_max = inf

    py_min = -inf
    py_max = inf

    pz_min = -inf
    pz_max = inf

    psi_min = -math.pi
    psi_max = math.pi

    v_max = 1.0
    w_max = 2.0
    a_max = 1.5
    aw_max = 8.0

    # Canonical start.
    px_cur = 0.0
    py_cur = 0.0
    pz_cur = 0.0
    psi_cur = 0.0

    vx_cur = 0.0
    vy_cur = 0.0
    vz_cur = 0.0
    w_cur = 0.0

    # Canonical target.
    pxf = 1.0
    pyf = 0.0
    pzf = 0.0
    psif = 0.0

    problem = lib.create_point_set_problem(
        N,
        tf,

        px_max,
        px_min,
        py_max,
        py_min,
        pz_max,
        pz_min,
        psi_max,
        psi_min,

        v_max,
        w_max,
        a_max,
        aw_max,

        px_cur,
        py_cur,
        pz_cur,
        psi_cur,

        vx_cur,
        vy_cur,
        vz_cur,
        w_cur,

        pxf,
        pyf,
        pzf,
        psif,

        sphere_x,
        sphere_y,
        sphere_z,
        sphere_radius,
    )

    if not problem:
        raise RuntimeError(
            "create_point_set_problem() returned nullptr."
        )

    try:
        solve_start = time.perf_counter()

        solve_ok = lib.solve_point_set_problem(
            problem
        )

        solve_seconds = (
            time.perf_counter()
            - solve_start
        )

        if not solve_ok:
            raise RuntimeError(
                "C++ IPOPT solve failed for this case."
            )

        solution_size = (
            lib.get_solution_size(problem)
        )

        expected_solution_size = (
            12 * (N + 1)
        )

        if (
            solution_size
            != expected_solution_size
        ):
            raise RuntimeError(
                "Unexpected NLP solution size: {} "
                "(expected {}).".format(
                    solution_size,
                    expected_solution_size,
                )
            )

        buffer_type = (
            ctypes.c_double
            * solution_size
        )

        solution_buffer = buffer_type()

        copied = lib.get_solution(
            problem,
            solution_buffer,
            solution_size,
        )

        if copied != solution_size:
            raise RuntimeError(
                "get_solution copied {} values, "
                "expected {}.".format(
                    copied,
                    solution_size,
                )
            )

        solution = np.ctypeslib.as_array(
            solution_buffer
        ).copy()

        objective = (
            lib.get_final_objective_value(
                problem
            )
        )

    finally:
        lib.destroy_point_set_problem(
            problem
        )

    L = N + 1

    nlp = {
        "px": solution[0 * L:1 * L],
        "py": solution[1 * L:2 * L],
        "pz": solution[2 * L:3 * L],
        "psi": solution[3 * L:4 * L],
    }

    return nlp, objective, solve_seconds


# ============================================================
# Load the trained NN checkpoint.
# ============================================================
def load_nn(checkpoint_path):
    checkpoint_path = Path(
        checkpoint_path
    ).resolve()

    if not checkpoint_path.exists():
        raise FileNotFoundError(
            "Could not find trained model: {}".format(
                checkpoint_path
            )
        )

    checkpoint = torch.load(
        checkpoint_path,
        map_location="cpu",
    )

    model = BernsteinMLP(
        input_dim=int(
            checkpoint["input_dim"]
        ),
        output_dim=int(
            checkpoint["output_dim"]
        ),
        width=int(
            checkpoint.get(
                "width",
                128,
            )
        ),
        residual_blocks=int(
            checkpoint.get(
                "residual_blocks",
                3,
            )
        ),
    )

    model.load_state_dict(
        checkpoint["model_state_dict"]
    )

    model.eval()

    return model, checkpoint


# ============================================================
# Predict with the NN.
#
# Training preprocessing:
#
#   rho   = sqrt(y_o^2 + z_o^2)
#   theta = atan2(z_o, y_o)
#
# Network input:
#
#   [x_o, rho, r]
#
# Network output is in obstacle-aligned coordinates:
#
#   [Px, Py_aligned, Pz_aligned, Psi]
#
# Then rotate predicted Py/Pz back by +theta.
# ============================================================
def predict_nn(
    model,
    checkpoint,
    sphere_x,
    sphere_y,
    sphere_z,
    sphere_radius,
):
    rho = math.sqrt(
        sphere_y * sphere_y
        + sphere_z * sphere_z
    )

    if rho < 1e-14:
        theta = 0.0
    else:
        theta = math.atan2(
            sphere_z,
            sphere_y,
        )

    X = np.array(
        [
            sphere_x,
            rho,
            sphere_radius,
        ],
        dtype=np.float64,
    )

    X_mean = np.asarray(
        checkpoint["X_mean"],
        dtype=np.float64,
    )

    X_std = np.asarray(
        checkpoint["X_std"],
        dtype=np.float64,
    )

    Y_mean = np.asarray(
        checkpoint["Y_mean"],
        dtype=np.float64,
    )

    Y_std = np.asarray(
        checkpoint["Y_std"],
        dtype=np.float64,
    )

    X_norm = (
        (X - X_mean)
        / X_std
    )

    x_tensor = torch.tensor(
        X_norm[None, :],
        dtype=torch.float32,
    )

    start = time.perf_counter()

    with torch.no_grad():
        y_norm = (
            model(x_tensor)
            .cpu()
            .numpy()[0]
        )

    inference_seconds = (
        time.perf_counter()
        - start
    )

    y_aligned = (
        y_norm * Y_std
        + Y_mean
    )

    N = int(checkpoint["N"])
    L = N + 1

    px = y_aligned[
        0 * L:
        1 * L
    ].copy()

    py_aligned = y_aligned[
        1 * L:
        2 * L
    ].copy()

    pz_aligned = y_aligned[
        2 * L:
        3 * L
    ].copy()

    psi = y_aligned[
        3 * L:
        4 * L
    ].copy()

    c = math.cos(theta)
    s = math.sin(theta)

    # Inverse of the training rotation.
    py = (
        c * py_aligned
        - s * pz_aligned
    )

    pz = (
        s * py_aligned
        + c * pz_aligned
    )

    prediction = {
        "px": px,
        "py": py,
        "pz": pz,
        "psi": psi,
        "py_aligned": py_aligned,
        "pz_aligned": pz_aligned,
    }

    geometry = {
        "rho": rho,
        "theta": theta,
    }

    return (
        prediction,
        geometry,
        inference_seconds,
    )


# ============================================================
# Bernstein curve reconstruction.
# ============================================================
def bernstein_curve(
    control_points,
    samples=400,
):
    cp = np.asarray(
        control_points,
        dtype=np.float64,
    )

    N = len(cp) - 1

    s = np.linspace(
        0.0,
        1.0,
        samples,
    )

    curve = np.zeros_like(s)

    for i in range(N + 1):
        curve += (
            math.comb(N, i)
            * (s ** i)
            * ((1.0 - s) ** (N - i))
            * cp[i]
        )

    return curve


# ============================================================
# Dense sphere clearance.
# ============================================================
def compute_clearance(
    control_points,
    sphere_x,
    sphere_y,
    sphere_z,
    sphere_radius,
):
    x = bernstein_curve(
        control_points["px"]
    )

    y = bernstein_curve(
        control_points["py"]
    )

    z = bernstein_curve(
        control_points["pz"]
    )

    distance = np.sqrt(
        (x - sphere_x) ** 2
        + (y - sphere_y) ** 2
        + (z - sphere_z) ** 2
    )

    min_distance = float(
        np.min(distance)
    )

    clearance = (
        min_distance
        - sphere_radius
    )

    return min_distance, clearance


# ============================================================
# Control-point comparison table.
# ============================================================
def build_comparison_table(
    nlp,
    nn_prediction,
):
    rows = []

    for variable in [
        "px",
        "py",
        "pz",
        "psi",
    ]:
        for i in range(
            len(nlp[variable])
        ):
            nlp_value = float(
                nlp[variable][i]
            )

            nn_value = float(
                nn_prediction[variable][i]
            )

            rows.append(
                {
                    "variable": variable,
                    "cp_index": i,
                    "nlp": nlp_value,
                    "nn": nn_value,
                    "error_nn_minus_nlp": (
                        nn_value - nlp_value
                    ),
                    "absolute_error": abs(
                        nn_value - nlp_value
                    ),
                }
            )

    return pd.DataFrame(rows)


# ============================================================
# Plot NLP and NN trajectories with sphere.
# ============================================================
def save_trajectory_plot(
    nlp,
    nn_prediction,
    sphere_x,
    sphere_y,
    sphere_z,
    sphere_radius,
    output_path,
):
    nlp_x = bernstein_curve(
        nlp["px"]
    )

    nlp_y = bernstein_curve(
        nlp["py"]
    )

    nlp_z = bernstein_curve(
        nlp["pz"]
    )

    nn_x = bernstein_curve(
        nn_prediction["px"]
    )

    nn_y = bernstein_curve(
        nn_prediction["py"]
    )

    nn_z = bernstein_curve(
        nn_prediction["pz"]
    )

    u = np.linspace(
        0.0,
        2.0 * np.pi,
        42,
    )

    v = np.linspace(
        0.0,
        np.pi,
        28,
    )

    surface_x = (
        sphere_x
        + sphere_radius
        * np.outer(
            np.cos(u),
            np.sin(v),
        )
    )

    surface_y = (
        sphere_y
        + sphere_radius
        * np.outer(
            np.sin(u),
            np.sin(v),
        )
    )

    surface_z = (
        sphere_z
        + sphere_radius
        * np.outer(
            np.ones_like(u),
            np.cos(v),
        )
    )

    fig = plt.figure(
        figsize=(10, 8)
    )

    ax = fig.add_subplot(
        111,
        projection="3d",
    )

    ax.plot(
        nlp_x,
        nlp_y,
        nlp_z,
        linewidth=2.8,
        label="C++ NLP",
    )

    ax.plot(
        nn_x,
        nn_y,
        nn_z,
        linewidth=2.4,
        linestyle="--",
        label="Neural network",
    )

    ax.plot(
        [0.0, 1.0],
        [0.0, 0.0],
        [0.0, 0.0],
        linestyle=":",
        linewidth=1.4,
        label="Direct start-target line",
    )

    ax.plot_wireframe(
        surface_x,
        surface_y,
        surface_z,
        rstride=2,
        cstride=2,
        linewidth=0.55,
        alpha=0.55,
    )

    ax.scatter(
        [0.0],
        [0.0],
        [0.0],
        s=130,
        marker="*",
        label="Start",
    )

    ax.scatter(
        [1.0],
        [0.0],
        [0.0],
        s=110,
        marker="X",
        label="Target",
    )

    ax.scatter(
        [sphere_x],
        [sphere_y],
        [sphere_z],
        s=65,
        label="Sphere center",
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
        "C++ NLP vs trained neural network"
    )

    ax.legend(
        fontsize=8
    )

    fig.tight_layout()

    fig.savefig(
        output_path,
        dpi=220,
        bbox_inches="tight",
    )

    plt.close(fig)


# ============================================================
# Main
# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare C++ IPOPT/BeBOT control points with "
            "the trained PyTorch network for exactly the same "
            "canonical one-sphere planning problem."
        )
    )

    parser.add_argument(
        "--library",
        default="./libsingle_drone_single_obs_nlp_fast.so",
        help="Fast C++ NLP shared library",
    )

    parser.add_argument(
        "--model",
        default="./single_obs_mlp_best.pt",
        help="Trained PyTorch checkpoint",
    )

    # Default is an unseen blocking example.
    parser.add_argument(
        "--sphere-x",
        type=float,
        default=0.55,
    )

    parser.add_argument(
        "--sphere-y",
        type=float,
        default=0.08,
    )

    parser.add_argument(
        "--sphere-z",
        type=float,
        default=-0.05,
    )

    parser.add_argument(
        "--radius",
        type=float,
        default=0.14,
    )

    args = parser.parse_args()

    sphere_x = args.sphere_x
    sphere_y = args.sphere_y
    sphere_z = args.sphere_z
    sphere_radius = args.radius

    rho = math.sqrt(
        sphere_y * sphere_y
        + sphere_z * sphere_z
    )

    line_gap = (
        rho - sphere_radius
    )

    print(
        "============================================================"
    )
    print(
        "C++ NLP vs TRAINED NN"
    )
    print(
        "============================================================"
    )
    print(
        "Canonical start : [0, 0, 0]"
    )
    print(
        "Canonical target: [1, 0, 0]"
    )
    print(
        "Sphere center   : [{:.8f}, {:.8f}, {:.8f}]".format(
            sphere_x,
            sphere_y,
            sphere_z,
        )
    )
    print(
        "Sphere radius   : {:.8f}".format(
            sphere_radius
        )
    )
    print(
        "rho             : {:.8f}".format(
            rho
        )
    )
    print(
        "rho - radius    : {:.8f}".format(
            line_gap
        )
    )

    if line_gap < 0.0:
        print(
            "Geometry class  : BLOCKING"
        )
    else:
        ratio = (
            rho / sphere_radius
            if sphere_radius > 0.0
            else float("inf")
        )

        if ratio <= 1.50:
            print(
                "Geometry class  : NEAR-MISS"
            )
        else:
            print(
                "Geometry class  : CLEAR"
            )

    print(
        "============================================================\n"
    )

    # --------------------------------------------------------
    # C++ NLP
    # --------------------------------------------------------
    lib = load_nlp_library(
        args.library
    )

    nlp, objective, nlp_seconds = (
        solve_cpp_nlp(
            lib,
            sphere_x,
            sphere_y,
            sphere_z,
            sphere_radius,
        )
    )

    # --------------------------------------------------------
    # Trained NN
    # --------------------------------------------------------
    model, checkpoint = load_nn(
        args.model
    )

    (
        nn_prediction,
        geometry,
        nn_seconds,
    ) = predict_nn(
        model,
        checkpoint,
        sphere_x,
        sphere_y,
        sphere_z,
        sphere_radius,
    )

    # --------------------------------------------------------
    # Compare control points.
    # --------------------------------------------------------
    table = build_comparison_table(
        nlp,
        nn_prediction,
    )

    print(
        "\nCONTROL-POINT COMPARISON"
    )
    print(
        "------------------------------------------------------------"
    )

    with pd.option_context(
        "display.max_rows",
        100,
        "display.width",
        160,
        "display.precision",
        10,
    ):
        print(
            table.to_string(
                index=False
            )
        )

    overall_rmse = float(
        np.sqrt(
            np.mean(
                table[
                    "error_nn_minus_nlp"
                ].to_numpy() ** 2
            )
        )
    )

    overall_mae = float(
        table[
            "absolute_error"
        ].mean()
    )

    print(
        "\nControl-point RMSE = {:.10e}".format(
            overall_rmse
        )
    )

    print(
        "Control-point MAE  = {:.10e}".format(
            overall_mae
        )
    )

    print(
        "\nPer-variable comparison:"
    )

    for variable in [
        "px",
        "py",
        "pz",
        "psi",
    ]:
        subset = table[
            table["variable"]
            == variable
        ]

        rmse = float(
            np.sqrt(
                np.mean(
                    subset[
                        "error_nn_minus_nlp"
                    ].to_numpy() ** 2
                )
            )
        )

        mae = float(
            subset[
                "absolute_error"
            ].mean()
        )

        print(
            "  {:4s}: RMSE = {:.10e}, MAE = {:.10e}".format(
                variable,
                rmse,
                mae,
            )
        )

    # --------------------------------------------------------
    # Geometric comparison.
    # --------------------------------------------------------
    nlp_min_distance, nlp_clearance = (
        compute_clearance(
            nlp,
            sphere_x,
            sphere_y,
            sphere_z,
            sphere_radius,
        )
    )

    nn_min_distance, nn_clearance = (
        compute_clearance(
            nn_prediction,
            sphere_x,
            sphere_y,
            sphere_z,
            sphere_radius,
        )
    )

    print(
        "\nGEOMETRIC COMPARISON"
    )
    print(
        "------------------------------------------------------------"
    )
    print(
        "NLP min distance to sphere center : {:.10f}".format(
            nlp_min_distance
        )
    )
    print(
        "NLP clearance                     : {:.10f}".format(
            nlp_clearance
        )
    )
    print(
        "NN min distance to sphere center  : {:.10f}".format(
            nn_min_distance
        )
    )
    print(
        "NN clearance                      : {:.10f}".format(
            nn_clearance
        )
    )

    print(
        "\nTIMING"
    )
    print(
        "------------------------------------------------------------"
    )
    print(
        "C++ NLP solve wall time : {:.6f} s".format(
            nlp_seconds
        )
    )
    print(
        "NN forward-pass time    : {:.9f} s".format(
            nn_seconds
        )
    )

    if nn_seconds > 0.0:
        print(
            "Approx. NLP / NN ratio  : {:.1f}x".format(
                nlp_seconds
                / nn_seconds
            )
        )

    print(
        "\nNLP objective = {:.10f}".format(
            objective
        )
    )

    # --------------------------------------------------------
    # Save files without overwriting earlier comparisons.
    # --------------------------------------------------------
    csv_path = non_overwriting_path(
        "nlp_vs_nn_control_points.csv"
    )

    plot_path = non_overwriting_path(
        "nlp_vs_nn_trajectory.png"
    )

    table.to_csv(
        csv_path,
        index=False,
    )

    save_trajectory_plot(
        nlp,
        nn_prediction,
        sphere_x,
        sphere_y,
        sphere_z,
        sphere_radius,
        plot_path,
    )

    print(
        "\nSaved comparison CSV : {}".format(
            csv_path
        )
    )
    print(
        "Saved trajectory plot: {}".format(
            plot_path
        )
    )

    print(
        "============================================================"
    )


if __name__ == "__main__":
    main()

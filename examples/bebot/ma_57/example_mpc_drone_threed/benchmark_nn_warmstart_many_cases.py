#!/usr/bin/env python3

import argparse
import ctypes
import math
import random
import statistics
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from torch import nn


# ============================================================
# Neural network: MUST match training architecture
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
# Shared library interface
# ============================================================
def load_library(path):
    path = Path(path).resolve()

    if not path.exists():
        raise FileNotFoundError(
            f"Could not find library: {path}"
        )

    lib = ctypes.CDLL(str(path))

    lib.create_point_set_problem.argtypes = (
        [ctypes.c_int]
        + [ctypes.c_double] * 29
    )
    lib.create_point_set_problem.restype = ctypes.c_void_p

    lib.solve_point_set_problem_quiet.argtypes = [
        ctypes.c_void_p
    ]
    lib.solve_point_set_problem_quiet.restype = ctypes.c_int

    lib.set_constant_initial_guess.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
    ]
    lib.set_constant_initial_guess.restype = ctypes.c_int

    lib.set_position_yaw_initial_guess.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
    ]
    lib.set_position_yaw_initial_guess.restype = ctypes.c_int

    lib.get_solution_size.argtypes = [
        ctypes.c_void_p
    ]
    lib.get_solution_size.restype = ctypes.c_int

    lib.get_solution.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
    ]
    lib.get_solution.restype = ctypes.c_int

    lib.get_final_objective_value.argtypes = [
        ctypes.c_void_p
    ]
    lib.get_final_objective_value.restype = ctypes.c_double

    lib.destroy_point_set_problem.argtypes = [
        ctypes.c_void_p
    ]
    lib.destroy_point_set_problem.restype = None

    return lib


# ============================================================
# Canonical OCP
# ============================================================
def create_problem(
    lib,
    sphere_x,
    sphere_y,
    sphere_z,
    sphere_radius,
):
    N = 4
    tf = 0.9

    inf = float("inf")

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

    # Canonical start
    px_cur = 0.0
    py_cur = 0.0
    pz_cur = 0.0
    psi_cur = 0.0

    vx_cur = 0.0
    vy_cur = 0.0
    vz_cur = 0.0
    w_cur = 0.0

    # Canonical target
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
            "create_point_set_problem returned nullptr."
        )

    return problem


# ============================================================
# NN loading / inference
# ============================================================
def load_model(path):
    path = Path(path).resolve()

    if not path.exists():
        raise FileNotFoundError(
            f"Could not find model: {path}"
        )

    checkpoint = torch.load(
        path,
        map_location="cpu",
        weights_only=False,
    )

    model = BernsteinMLP(
        input_dim=int(checkpoint["input_dim"]),
        output_dim=int(checkpoint["output_dim"]),
        width=int(checkpoint.get("width", 128)),
        residual_blocks=int(
            checkpoint.get("residual_blocks", 3)
        ),
    )

    model.load_state_dict(
        checkpoint["model_state_dict"]
    )
    model.eval()

    return model, checkpoint


def predict_nn_guess(
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

    theta = (
        0.0
        if rho < 1e-14
        else math.atan2(
            sphere_z,
            sphere_y,
        )
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

    Xn = (
        (X - X_mean)
        / X_std
    )

    xt = torch.tensor(
        Xn[None, :],
        dtype=torch.float32,
    )

    start = time.perf_counter()

    with torch.no_grad():
        yn = (
            model(xt)
            .cpu()
            .numpy()[0]
        )

    inference_time = (
        time.perf_counter()
        - start
    )

    aligned = (
        yn * Y_std
        + Y_mean
    )

    N = int(checkpoint["N"])
    L = N + 1

    px = aligned[
        0 * L:
        1 * L
    ].copy()

    py_aligned = aligned[
        1 * L:
        2 * L
    ].copy()

    pz_aligned = aligned[
        2 * L:
        3 * L
    ].copy()

    psi = aligned[
        3 * L:
        4 * L
    ].copy()

    c = math.cos(theta)
    s = math.sin(theta)

    py = (
        c * py_aligned
        - s * pz_aligned
    )

    pz = (
        s * py_aligned
        + c * pz_aligned
    )

    position_yaw = np.concatenate(
        [
            px,
            py,
            pz,
            psi,
        ]
    ).astype(
        np.float64,
        copy=False,
    )

    return (
        position_yaw,
        inference_time,
    )


# ============================================================
# Solution helpers
# ============================================================
def get_solution(lib, problem):
    n = lib.get_solution_size(
        problem
    )

    if n <= 0:
        return None, float("nan")

    buffer_type = (
        ctypes.c_double
        * n
    )

    buffer = buffer_type()

    copied = lib.get_solution(
        problem,
        buffer,
        n,
    )

    if copied != n:
        return None, float("nan")

    solution = np.ctypeslib.as_array(
        buffer
    ).copy()

    objective = float(
        lib.get_final_objective_value(
            problem
        )
    )

    return solution, objective


def solve_standard(
    lib,
    sphere,
    constant_guess,
):
    problem = create_problem(
        lib,
        *sphere,
    )

    try:
        ok = lib.set_constant_initial_guess(
            problem,
            float(constant_guess),
        )

        if not ok:
            return None

        start = time.perf_counter()

        solve_ok = (
            lib.solve_point_set_problem_quiet(
                problem
            )
        )

        solve_time = (
            time.perf_counter()
            - start
        )

        if not solve_ok:
            return {
                "success": False,
                "solve_time": solve_time,
                "objective": float("nan"),
                "solution": None,
            }

        solution, objective = get_solution(
            lib,
            problem,
        )

        return {
            "success": solution is not None,
            "solve_time": solve_time,
            "objective": objective,
            "solution": solution,
        }

    finally:
        lib.destroy_point_set_problem(
            problem
        )


def solve_nn_warm(
    lib,
    sphere,
    nn_guess,
):
    problem = create_problem(
        lib,
        *sphere,
    )

    try:
        guess = np.ascontiguousarray(
            nn_guess,
            dtype=np.float64,
        )

        ok = (
            lib.set_position_yaw_initial_guess(
                problem,
                guess.ctypes.data_as(
                    ctypes.POINTER(
                        ctypes.c_double
                    )
                ),
                int(guess.size),
            )
        )

        if not ok:
            return None

        start = time.perf_counter()

        solve_ok = (
            lib.solve_point_set_problem_quiet(
                problem
            )
        )

        solve_time = (
            time.perf_counter()
            - start
        )

        if not solve_ok:
            return {
                "success": False,
                "solve_time": solve_time,
                "objective": float("nan"),
                "solution": None,
            }

        solution, objective = get_solution(
            lib,
            problem,
        )

        return {
            "success": solution is not None,
            "solve_time": solve_time,
            "objective": objective,
            "solution": solution,
        }

    finally:
        lib.destroy_point_set_problem(
            problem
        )


# ============================================================
# Geometry / clearance
# ============================================================
def bernstein_curve(
    cp,
    samples=300,
):
    cp = np.asarray(
        cp,
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


def clearance(
    solution,
    sphere,
):
    if solution is None:
        return float("nan")

    sx, sy, sz, r = sphere

    L = 5

    px = solution[
        0 * L:
        1 * L
    ]
    py = solution[
        1 * L:
        2 * L
    ]
    pz = solution[
        2 * L:
        3 * L
    ]

    x = bernstein_curve(px)
    y = bernstein_curve(py)
    z = bernstein_curve(pz)

    d = np.sqrt(
        (x - sx) ** 2
        + (y - sy) ** 2
        + (z - sz) ** 2
    )

    return float(
        d.min()
        - r
    )


# ============================================================
# Sampling: same geometry distribution used for training
# ============================================================
def sample_case(
    rng,
    case_type,
):
    radius = rng.uniform(
        0.06,
        0.20,
    )

    margin = 0.05

    x_low = (
        radius + margin
    )

    x_high = (
        1.0
        - radius
        - margin
    )

    sphere_x = rng.uniform(
        x_low,
        x_high,
    )

    if case_type == "blocking":
        rho = rng.uniform(
            0.05 * radius,
            0.90 * radius,
        )

    elif case_type == "near_miss":
        rho = rng.uniform(
            1.05 * radius,
            1.50 * radius,
        )

    elif case_type == "clear":
        rho = rng.uniform(
            1.75 * radius,
            3.00 * radius,
        )

    else:
        raise ValueError(
            f"Unknown case type: {case_type}"
        )

    theta = rng.uniform(
        0.0,
        2.0 * math.pi,
    )

    sphere_y = (
        rho
        * math.cos(theta)
    )

    sphere_z = (
        rho
        * math.sin(theta)
    )

    return (
        sphere_x,
        sphere_y,
        sphere_z,
        radius,
        rho,
        theta,
    )


def build_case_schedule(
    n_cases,
    rng,
):
    n_blocking = int(
        round(
            0.50 * n_cases
        )
    )

    n_near = int(
        round(
            0.25 * n_cases
        )
    )

    n_clear = (
        n_cases
        - n_blocking
        - n_near
    )

    schedule = (
        ["blocking"] * n_blocking
        + ["near_miss"] * n_near
        + ["clear"] * n_clear
    )

    rng.shuffle(
        schedule
    )

    return schedule


# ============================================================
# Statistics / plots
# ============================================================
def safe_median(values):
    values = [
        float(v)
        for v in values
        if np.isfinite(v)
    ]

    if not values:
        return float("nan")

    return statistics.median(
        values
    )


def safe_mean(values):
    values = [
        float(v)
        for v in values
        if np.isfinite(v)
    ]

    if not values:
        return float("nan")

    return statistics.mean(
        values
    )


def print_category_summary(
    df,
    category,
):
    if category == "ALL":
        sub = df.copy()
    else:
        sub = df[
            df["case_type"]
            == category
        ].copy()

    matched = sub[
        sub["standard_success"]
        & sub["warm_success"]
    ].copy()

    if len(matched) == 0:
        print(
            f"{category:10s}: no matched successful cases"
        )
        return

    standard_ms = (
        matched[
            "standard_solve_ms"
        ].to_numpy()
    )

    warm_ms = (
        matched[
            "warm_solve_ms"
        ].to_numpy()
    )

    end_to_end_ms = (
        matched[
            "warm_end_to_end_ms"
        ].to_numpy()
    )

    speedup = (
        standard_ms
        / end_to_end_ms
    )

    warm_faster_fraction = float(
        np.mean(
            end_to_end_ms
            < standard_ms
        )
    )

    print(
        (
            f"{category:10s}: "
            f"matched={len(matched):4d} | "
            f"std med={np.median(standard_ms):8.4f} ms | "
            f"NN+NLP med={np.median(end_to_end_ms):8.4f} ms | "
            f"median speedup={np.median(speedup):6.3f}x | "
            f"warm faster={100.0 * warm_faster_fraction:6.2f}%"
        )
    )


def save_timing_scatter(
    df,
    output_path,
):
    matched = df[
        df["standard_success"]
        & df["warm_success"]
    ].copy()

    fig, ax = plt.subplots(
        figsize=(8, 7)
    )

    for case_type in [
        "blocking",
        "near_miss",
        "clear",
    ]:
        sub = matched[
            matched["case_type"]
            == case_type
        ]

        if len(sub) == 0:
            continue

        ax.scatter(
            sub["standard_solve_ms"],
            sub["warm_end_to_end_ms"],
            s=28,
            alpha=0.65,
            label=case_type,
        )

    if len(matched) > 0:
        low = float(
            min(
                matched[
                    "standard_solve_ms"
                ].min(),
                matched[
                    "warm_end_to_end_ms"
                ].min(),
            )
        )

        high = float(
            max(
                matched[
                    "standard_solve_ms"
                ].max(),
                matched[
                    "warm_end_to_end_ms"
                ].max(),
            )
        )

        ax.plot(
            [low, high],
            [low, high],
            linestyle="--",
            linewidth=1.4,
            label="Equal time",
        )

    ax.set_xlabel(
        "Standard NLP solve time [ms]"
    )

    ax.set_ylabel(
        "NN inference + warm-start NLP [ms]"
    )

    ax.set_title(
        "Per-case computation time"
    )

    ax.grid(
        alpha=0.25
    )

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        output_path,
        dpi=220,
        bbox_inches="tight",
    )

    plt.close(fig)


def save_speedup_histogram(
    df,
    output_path,
):
    matched = df[
        df["standard_success"]
        & df["warm_success"]
    ].copy()

    if len(matched) == 0:
        return

    speedup = (
        matched[
            "standard_solve_ms"
        ].to_numpy()
        / matched[
            "warm_end_to_end_ms"
        ].to_numpy()
    )

    fig, ax = plt.subplots(
        figsize=(8, 6)
    )

    ax.hist(
        speedup,
        bins=30,
        alpha=0.8,
    )

    ax.axvline(
        1.0,
        linestyle="--",
        linewidth=1.5,
        label="No speedup",
    )

    ax.set_xlabel(
        "End-to-end speedup factor"
    )

    ax.set_ylabel(
        "Number of cases"
    )

    ax.set_title(
        "NN warm-start speedup distribution"
    )

    ax.grid(
        alpha=0.25
    )

    ax.legend()

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
            "Benchmark standard NLP initialization versus "
            "NN-informed NLP initialization over many random "
            "single-sphere obstacle cases."
        )
    )

    parser.add_argument(
        "--library",
        default="./libsingle_drone_single_obs_nlp_warmstart.so",
    )

    parser.add_argument(
        "--model",
        default="./single_obs_mlp_best.pt",
    )

    parser.add_argument(
        "--cases",
        type=int,
        default=100,
        help="Number of random obstacle cases.",
    )

    parser.add_argument(
        "--standard-guess",
        choices=[
            "zero",
            "one",
        ],
        default="zero",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=20260821,
    )

    parser.add_argument(
        "--output-prefix",
        default="nn_warmstart_many_cases",
    )

    args = parser.parse_args()

    if args.cases <= 0:
        raise ValueError(
            "--cases must be positive."
        )

    # More stable timing for a tiny CPU model.
    torch.set_num_threads(1)

    rng = random.Random(
        args.seed
    )

    lib = load_library(
        args.library
    )

    model, checkpoint = load_model(
        args.model
    )

    # Warm up PyTorch once before real timing.
    warmup_case = (
        0.55,
        0.08,
        -0.05,
        0.14,
    )

    for _ in range(50):
        predict_nn_guess(
            model,
            checkpoint,
            *warmup_case,
        )

    schedule = build_case_schedule(
        args.cases,
        rng,
    )

    constant_guess = (
        0.0
        if args.standard_guess
        == "zero"
        else 1.0
    )

    rows = []

    print(
        "============================================================"
    )
    print(
        "MANY-CASE NN WARM-START BENCHMARK"
    )
    print(
        "============================================================"
    )
    print(
        f"Cases          : {args.cases}"
    )
    print(
        f"Standard guess : all {args.standard_guess}"
    )
    print(
        f"Seed           : {args.seed}"
    )
    print(
        "Distribution   : 50% blocking, "
        "25% near_miss, 25% clear"
    )
    print(
        "============================================================"
    )

    for case_id, case_type in enumerate(
        schedule,
        start=1,
    ):
        (
            sphere_x,
            sphere_y,
            sphere_z,
            radius,
            rho,
            theta,
        ) = sample_case(
            rng,
            case_type,
        )

        sphere = (
            sphere_x,
            sphere_y,
            sphere_z,
            radius,
        )

        # ----------------------------------------------------
        # NN prediction
        # ----------------------------------------------------
        nn_guess, nn_time = (
            predict_nn_guess(
                model,
                checkpoint,
                *sphere,
            )
        )

        # ----------------------------------------------------
        # Standard NLP
        # ----------------------------------------------------
        standard = solve_standard(
            lib,
            sphere,
            constant_guess,
        )

        # ----------------------------------------------------
        # NN-initialized NLP
        # ----------------------------------------------------
        warm = solve_nn_warm(
            lib,
            sphere,
            nn_guess,
        )

        standard_success = (
            standard is not None
            and standard["success"]
        )

        warm_success = (
            warm is not None
            and warm["success"]
        )

        standard_time_ms = (
            standard["solve_time"] * 1e3
            if standard is not None
            else float("nan")
        )

        warm_time_ms = (
            warm["solve_time"] * 1e3
            if warm is not None
            else float("nan")
        )

        nn_time_ms = (
            nn_time * 1e3
        )

        warm_end_to_end_ms = (
            warm_time_ms
            + nn_time_ms
            if warm is not None
            else float("nan")
        )

        standard_objective = (
            standard["objective"]
            if standard_success
            else float("nan")
        )

        warm_objective = (
            warm["objective"]
            if warm_success
            else float("nan")
        )

        standard_clearance = (
            clearance(
                standard["solution"],
                sphere,
            )
            if standard_success
            else float("nan")
        )

        warm_clearance = (
            clearance(
                warm["solution"],
                sphere,
            )
            if warm_success
            else float("nan")
        )

        if (
            standard_success
            and warm_success
        ):
            speedup = (
                standard_time_ms
                / warm_end_to_end_ms
            )

            objective_relative_difference = (
                abs(
                    warm_objective
                    - standard_objective
                )
                / max(
                    abs(
                        standard_objective
                    ),
                    1e-12,
                )
            )
        else:
            speedup = float("nan")
            objective_relative_difference = (
                float("nan")
            )

        rows.append(
            {
                "case_id": case_id,
                "case_type": case_type,

                "sphere_x": sphere_x,
                "sphere_y": sphere_y,
                "sphere_z": sphere_z,
                "sphere_radius": radius,
                "rho": rho,
                "rho_over_r": (
                    rho / radius
                ),
                "theta": theta,

                "standard_guess": (
                    args.standard_guess
                ),

                "standard_success": (
                    standard_success
                ),
                "warm_success": (
                    warm_success
                ),

                "nn_inference_ms": (
                    nn_time_ms
                ),

                "standard_solve_ms": (
                    standard_time_ms
                ),

                "warm_solve_ms": (
                    warm_time_ms
                ),

                "warm_end_to_end_ms": (
                    warm_end_to_end_ms
                ),

                "end_to_end_speedup": (
                    speedup
                ),

                "standard_objective": (
                    standard_objective
                ),

                "warm_objective": (
                    warm_objective
                ),

                "objective_relative_difference": (
                    objective_relative_difference
                ),

                "standard_clearance": (
                    standard_clearance
                ),

                "warm_clearance": (
                    warm_clearance
                ),
            }
        )

        if (
            case_id == 1
            or case_id % 10 == 0
            or case_id == args.cases
        ):
            print(
                f"Completed {case_id:4d} / {args.cases}"
            )

    df = pd.DataFrame(
        rows
    )

    csv_path = Path(
        args.output_prefix
        + ".csv"
    )

    scatter_path = Path(
        args.output_prefix
        + "_timing_scatter.png"
    )

    hist_path = Path(
        args.output_prefix
        + "_speedup_histogram.png"
    )

    df.to_csv(
        csv_path,
        index=False,
    )

    save_timing_scatter(
        df,
        scatter_path,
    )

    save_speedup_histogram(
        df,
        hist_path,
    )

    # ========================================================
    # Final summary
    # ========================================================
    n_standard_success = int(
        df[
            "standard_success"
        ].sum()
    )

    n_warm_success = int(
        df[
            "warm_success"
        ].sum()
    )

    both_success = df[
        df["standard_success"]
        & df["warm_success"]
    ].copy()

    standard_only = int(
        (
            df["standard_success"]
            & ~df["warm_success"]
        ).sum()
    )

    warm_only = int(
        (
            ~df["standard_success"]
            & df["warm_success"]
        ).sum()
    )

    print(
        "\n============================================================"
    )
    print(
        "FINAL COMPARISON"
    )
    print(
        "============================================================"
    )

    print(
        f"Standard NLP successes : "
        f"{n_standard_success} / {args.cases}"
    )

    print(
        f"NN-warm NLP successes  : "
        f"{n_warm_success} / {args.cases}"
    )

    print(
        f"Both succeeded         : "
        f"{len(both_success)} / {args.cases}"
    )

    print(
        f"Standard only success  : "
        f"{standard_only}"
    )

    print(
        f"NN-warm only success   : "
        f"{warm_only}"
    )

    print(
        "------------------------------------------------------------"
    )

    print_category_summary(
        df,
        "ALL",
    )

    print_category_summary(
        df,
        "blocking",
    )

    print_category_summary(
        df,
        "near_miss",
    )

    print_category_summary(
        df,
        "clear",
    )

    if len(both_success) > 0:
        speedups = (
            both_success[
                "end_to_end_speedup"
            ].to_numpy()
        )

        objective_diff = (
            both_success[
                "objective_relative_difference"
            ].to_numpy()
        )

        clearance_diff = np.abs(
            both_success[
                "warm_clearance"
            ].to_numpy()
            - both_success[
                "standard_clearance"
            ].to_numpy()
        )

        print(
            "------------------------------------------------------------"
        )

        print(
            "Median NN inference         : "
            f"{np.median(both_success['nn_inference_ms']):.6f} ms"
        )

        print(
            "Median standard NLP         : "
            f"{np.median(both_success['standard_solve_ms']):.6f} ms"
        )

        print(
            "Median NN-warm NLP only     : "
            f"{np.median(both_success['warm_solve_ms']):.6f} ms"
        )

        print(
            "Median NN + warm NLP        : "
            f"{np.median(both_success['warm_end_to_end_ms']):.6f} ms"
        )

        print(
            "Median end-to-end speedup   : "
            f"{np.median(speedups):.4f}x"
        )

        print(
            "Mean end-to-end speedup     : "
            f"{np.mean(speedups):.4f}x"
        )

        print(
            "Warm faster than standard   : "
            f"{100.0 * np.mean(speedups > 1.0):.2f}%"
        )

        print(
            "Median relative objective Δ : "
            f"{100.0 * np.median(objective_diff):.4f}%"
        )

        print(
            "Median |clearance Δ|        : "
            f"{np.median(clearance_diff):.8f}"
        )

        print(
            "Worst relative objective Δ  : "
            f"{100.0 * np.max(objective_diff):.4f}%"
        )

    print(
        "------------------------------------------------------------"
    )

    print(
        f"Saved CSV          : {csv_path}"
    )

    print(
        f"Saved timing plot  : {scatter_path}"
    )

    print(
        f"Saved speedup plot : {hist_path}"
    )

    print(
        "============================================================"
    )


if __name__ == "__main__":
    main()
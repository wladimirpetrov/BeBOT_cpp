#!/usr/bin/env python3

import argparse
import ctypes
import math
import statistics
import time
from pathlib import Path

import numpy as np
import torch
from torch import nn


# ============================================================
# Same network architecture used during training
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
# Shared-library interface
# ============================================================
def load_library(path):
    path = Path(path).resolve()

    if not path.exists():
        raise FileNotFoundError(
            "Could not find warm-start library: {}".format(path)
        )

    lib = ctypes.CDLL(str(path))

    lib.create_point_set_problem.argtypes = (
        [ctypes.c_int]
        + [ctypes.c_double] * 29
    )
    lib.create_point_set_problem.restype = ctypes.c_void_p

    lib.solve_point_set_problem.argtypes = [
        ctypes.c_void_p
    ]
    lib.solve_point_set_problem.restype = ctypes.c_int

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
# Canonical C++ OCP creation
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
            "create_point_set_problem() returned nullptr."
        )

    return problem


# ============================================================
# NN loading and prediction
# ============================================================
def load_model(path):
    path = Path(path).resolve()

    if not path.exists():
        raise FileNotFoundError(
            "Could not find model: {}".format(path)
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


def prepare_nn_input(
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

    X_normalized = (
        (X - X_mean)
        / X_std
    )

    tensor = torch.tensor(
        X_normalized[None, :],
        dtype=torch.float32,
    )

    return tensor, theta


def decode_nn_output(
    normalized_output,
    checkpoint,
    theta,
):
    Y_mean = np.asarray(
        checkpoint["Y_mean"],
        dtype=np.float64,
    )

    Y_std = np.asarray(
        checkpoint["Y_std"],
        dtype=np.float64,
    )

    aligned = (
        normalized_output
        * Y_std
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

    # 4*(N+1) array required by the C++ warm-start interface.
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

    return position_yaw


def predict_once(
    model,
    checkpoint,
    tensor,
    theta,
):
    with torch.no_grad():
        output = (
            model(tensor)
            .cpu()
            .numpy()[0]
        )

    return decode_nn_output(
        output,
        checkpoint,
        theta,
    )


# ============================================================
# C++ solution helpers
# ============================================================
def get_solution(lib, problem):
    n = lib.get_solution_size(
        problem
    )

    if n <= 0:
        raise RuntimeError(
            "No solution stored by NLP."
        )

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
        raise RuntimeError(
            "Failed to copy complete NLP solution."
        )

    solution = np.ctypeslib.as_array(
        buffer
    ).copy()

    objective = (
        lib.get_final_objective_value(
            problem
        )
    )

    return solution, objective


def solve_with_constant_guess(
    lib,
    sphere,
    guess_value,
    quiet=True,
):
    t_total_start = time.perf_counter()

    problem = create_problem(
        lib,
        *sphere,
    )

    set_ok = lib.set_constant_initial_guess(
        problem,
        float(guess_value),
    )

    if not set_ok:
        lib.destroy_point_set_problem(problem)
        raise RuntimeError(
            "Could not set constant initial guess."
        )

    t_solve_start = time.perf_counter()

    solve_ok = (
        lib.solve_point_set_problem_quiet(problem)
        if quiet
        else lib.solve_point_set_problem(problem)
    )

    solve_time = (
        time.perf_counter()
        - t_solve_start
    )

    if not solve_ok:
        lib.destroy_point_set_problem(problem)
        return None

    solution, objective = get_solution(
        lib,
        problem,
    )

    total_time = (
        time.perf_counter()
        - t_total_start
    )

    lib.destroy_point_set_problem(
        problem
    )

    return {
        "solution": solution,
        "objective": objective,
        "solve_time": solve_time,
        "total_time": total_time,
    }


def solve_with_nn_guess(
    lib,
    sphere,
    position_yaw_guess,
    quiet=True,
):
    t_total_start = time.perf_counter()

    problem = create_problem(
        lib,
        *sphere,
    )

    guess = np.ascontiguousarray(
        position_yaw_guess,
        dtype=np.float64,
    )

    set_ok = (
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

    if not set_ok:
        lib.destroy_point_set_problem(problem)

        raise RuntimeError(
            "Could not set NN position/yaw initial guess."
        )

    t_solve_start = time.perf_counter()

    solve_ok = (
        lib.solve_point_set_problem_quiet(problem)
        if quiet
        else lib.solve_point_set_problem(problem)
    )

    solve_time = (
        time.perf_counter()
        - t_solve_start
    )

    if not solve_ok:
        lib.destroy_point_set_problem(problem)
        return None

    solution, objective = get_solution(
        lib,
        problem,
    )

    total_time = (
        time.perf_counter()
        - t_total_start
    )

    lib.destroy_point_set_problem(
        problem
    )

    return {
        "solution": solution,
        "objective": objective,
        "solve_time": solve_time,
        "total_time": total_time,
    }


# ============================================================
# Geometry diagnostics
# ============================================================
def bernstein_curve(cp, samples=500):
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


def solution_clearance(
    solution,
    sphere_x,
    sphere_y,
    sphere_z,
    sphere_radius,
):
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

    distance = np.sqrt(
        (x - sphere_x) ** 2
        + (y - sphere_y) ** 2
        + (z - sphere_z) ** 2
    )

    return (
        float(distance.min())
        - sphere_radius
    )


# ============================================================
# Timing statistics
# ============================================================
def describe(values):
    values = list(values)

    return {
        "mean": statistics.mean(values),
        "median": statistics.median(values),
        "min": min(values),
        "max": max(values),
    }


def print_stats(label, values):
    stats = describe(values)

    print(
        "{:<32s} mean {:10.6f} ms | "
        "median {:10.6f} ms | "
        "min {:10.6f} ms | "
        "max {:10.6f} ms".format(
            label,
            stats["mean"] * 1e3,
            stats["median"] * 1e3,
            stats["min"] * 1e3,
            stats["max"] * 1e3,
        )
    )


# ============================================================
# Main benchmark
# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare IPOPT from a generic constant initial guess "
            "against IPOPT initialized by the trained neural network."
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

    parser.add_argument(
        "--standard-guess",
        choices=[
            "zero",
            "one",
        ],
        default="zero",
    )

    parser.add_argument(
        "--nlp-repeats",
        type=int,
        default=20,
    )

    parser.add_argument(
        "--nn-repeats",
        type=int,
        default=1000,
    )

    args = parser.parse_args()

    # Tiny MLPs often benchmark more cleanly with one CPU thread.
    torch.set_num_threads(1)

    lib = load_library(
        args.library
    )

    model, checkpoint = load_model(
        args.model
    )

    sphere = (
        args.sphere_x,
        args.sphere_y,
        args.sphere_z,
        args.radius,
    )

    tensor, theta = prepare_nn_input(
        checkpoint,
        args.sphere_x,
        args.sphere_y,
        args.sphere_z,
        args.radius,
    )

    # --------------------------------------------------------
    # NN warm-up: do not benchmark first-call initialization.
    # --------------------------------------------------------
    for _ in range(50):
        predict_once(
            model,
            checkpoint,
            tensor,
            theta,
        )

    nn_times = []
    nn_guess = None

    for _ in range(args.nn_repeats):
        start = time.perf_counter()

        nn_guess = predict_once(
            model,
            checkpoint,
            tensor,
            theta,
        )

        nn_times.append(
            time.perf_counter()
            - start
        )

    # --------------------------------------------------------
    # Benchmark generic constant initial guess.
    # --------------------------------------------------------
    constant_value = (
        0.0
        if args.standard_guess == "zero"
        else 1.0
    )

    standard_results = []

    for _ in range(args.nlp_repeats):
        result = solve_with_constant_guess(
            lib,
            sphere,
            constant_value,
            quiet=True,
        )

        if result is not None:
            standard_results.append(
                result
            )

    # --------------------------------------------------------
    # Benchmark NN-informed initial guess.
    # --------------------------------------------------------
    warm_results = []

    for _ in range(args.nlp_repeats):
        result = solve_with_nn_guess(
            lib,
            sphere,
            nn_guess,
            quiet=True,
        )

        if result is not None:
            warm_results.append(
                result
            )

    if not standard_results:
        raise RuntimeError(
            "All standard-initialization NLP solves failed."
        )

    if not warm_results:
        raise RuntimeError(
            "All NN-initialized NLP solves failed."
        )

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------
    rho = math.sqrt(
        args.sphere_y ** 2
        + args.sphere_z ** 2
    )

    print(
        "\n============================================================"
    )
    print(
        "NN-INFORMED NLP INITIALIZATION BENCHMARK"
    )
    print(
        "============================================================"
    )
    print(
        "Start              : [0, 0, 0]"
    )
    print(
        "Target             : [1, 0, 0]"
    )
    print(
        "Sphere             : "
        "[{:.6f}, {:.6f}, {:.6f}], r={:.6f}".format(
            args.sphere_x,
            args.sphere_y,
            args.sphere_z,
            args.radius,
        )
    )
    print(
        "rho/r              : {:.6f}".format(
            rho / args.radius
        )
    )
    print(
        "Standard guess     : all {}".format(
            args.standard_guess
        )
    )
    print(
        "Successful standard: {} / {}".format(
            len(standard_results),
            args.nlp_repeats,
        )
    )
    print(
        "Successful NN-warm : {} / {}".format(
            len(warm_results),
            args.nlp_repeats,
        )
    )
    print(
        "------------------------------------------------------------"
    )

    standard_solve_times = [
        r["solve_time"]
        for r in standard_results
    ]

    warm_solve_times = [
        r["solve_time"]
        for r in warm_results
    ]

    standard_total_times = [
        r["total_time"]
        for r in standard_results
    ]

    warm_total_times = [
        r["total_time"]
        for r in warm_results
    ]

    print_stats(
        "NN inference",
        nn_times,
    )

    print_stats(
        "NLP solve: standard guess",
        standard_solve_times,
    )

    print_stats(
        "NLP solve: NN initial guess",
        warm_solve_times,
    )

    print_stats(
        "C++ total: standard guess",
        standard_total_times,
    )

    print_stats(
        "C++ total: NN initial guess",
        warm_total_times,
    )

    standard_median = statistics.median(
        standard_solve_times
    )

    warm_median = statistics.median(
        warm_solve_times
    )

    nn_median = statistics.median(
        nn_times
    )

    print(
        "------------------------------------------------------------"
    )

    print(
        "Median NLP speedup from NN guess : {:.3f}x".format(
            standard_median
            / warm_median
        )
    )

    combined_warm_median = (
        nn_median
        + warm_median
    )

    print(
        "Median NN + warm NLP time        : {:.6f} ms".format(
            combined_warm_median
            * 1e3
        )
    )

    print(
        "Standard NLP median time         : {:.6f} ms".format(
            standard_median
            * 1e3
        )
    )

    print(
        "End-to-end median speedup        : {:.3f}x".format(
            standard_median
            / combined_warm_median
        )
    )

    # --------------------------------------------------------
    # Compare converged solution quality.
    # Use the last successful solve of each type.
    # --------------------------------------------------------
    standard = standard_results[-1]
    warm = warm_results[-1]

    L = 5

    standard_posyaw = np.concatenate(
        [
            standard["solution"][0 * L:1 * L],
            standard["solution"][1 * L:2 * L],
            standard["solution"][2 * L:3 * L],
            standard["solution"][3 * L:4 * L],
        ]
    )

    warm_posyaw = np.concatenate(
        [
            warm["solution"][0 * L:1 * L],
            warm["solution"][1 * L:2 * L],
            warm["solution"][2 * L:3 * L],
            warm["solution"][3 * L:4 * L],
        ]
    )

    cp_rmse = float(
        np.sqrt(
            np.mean(
                (
                    warm_posyaw
                    - standard_posyaw
                ) ** 2
            )
        )
    )

    standard_clearance = solution_clearance(
        standard["solution"],
        *sphere,
    )

    warm_clearance = solution_clearance(
        warm["solution"],
        *sphere,
    )

    print(
        "\nCONVERGED-SOLUTION CHECK"
    )
    print(
        "------------------------------------------------------------"
    )
    print(
        "Standard objective   : {:.12f}".format(
            standard["objective"]
        )
    )
    print(
        "NN-warm objective    : {:.12f}".format(
            warm["objective"]
        )
    )
    print(
        "Px/Py/Pz/Psi RMSE    : {:.12e}".format(
            cp_rmse
        )
    )
    print(
        "Standard clearance   : {:.12f}".format(
            standard_clearance
        )
    )
    print(
        "NN-warm clearance    : {:.12f}".format(
            warm_clearance
        )
    )
    print(
        "============================================================"
    )

    # --------------------------------------------------------
    # One visible IPOPT run for each initialization.
    # This is intentionally AFTER timing so console output does
    # not contaminate the benchmark.
    # --------------------------------------------------------
    print(
        "\nNow running ONE visible standard-guess solve "
        "and ONE visible NN-initialized solve so you can inspect "
        "the IPOPT iteration summaries.\n"
    )

    print(
        "================ STANDARD INITIAL GUESS ================"
    )
    visible_standard = solve_with_constant_guess(
        lib,
        sphere,
        constant_value,
        quiet=False,
    )

    print(
        "\n================ NN INITIAL GUESS ======================="
    )
    visible_warm = solve_with_nn_guess(
        lib,
        sphere,
        nn_guess,
        quiet=False,
    )


if __name__ == "__main__":
    main()

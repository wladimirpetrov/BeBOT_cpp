#!/usr/bin/env python3

import argparse
import json
import math
import random
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset


# ============================================================
# Reproducibility
# ============================================================
def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


# ============================================================
# Geometry normalization
#
# Dataset geometry:
#   start  = [0, 0, 0]
#   target = [1, 0, 0]
#
# Original obstacle:
#   [sphere_x, sphere_y, sphere_z], radius
#
# Let
#   rho   = sqrt(sphere_y^2 + sphere_z^2)
#   theta = atan2(sphere_z, sphere_y)
#
# Rotate the whole problem about the x-axis by -theta so that
# the obstacle is always:
#
#   [sphere_x, rho, 0]
#
# Network input therefore becomes:
#
#   [sphere_x, rho, sphere_radius]
#
# This removes an irrelevant rotational degree of freedom.
# ============================================================
def rotate_yz_to_obstacle_aligned(y, z, theta):
    c = np.cos(theta)
    s = np.sin(theta)

    y_aligned = c * y + s * z
    z_aligned = -s * y + c * z

    return y_aligned, z_aligned


def rotate_yz_from_obstacle_aligned(y_aligned, z_aligned, theta):
    c = np.cos(theta)
    s = np.sin(theta)

    y = c * y_aligned - s * z_aligned
    z = s * y_aligned + c * z_aligned

    return y, z


# ============================================================
# Dataset preparation
# ============================================================
def prepare_dataset(csv_path):
    df = pd.read_csv(csv_path)

    required = [
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
    ]

    for col in required:
        if col not in df.columns:
            raise ValueError(
                "Dataset is missing required column: {}".format(col)
            )

    unique_N = sorted(df["N"].unique().tolist())

    if len(unique_N) != 1:
        raise ValueError(
            "This script expects one Bernstein degree N. Found: {}".format(
                unique_N
            )
        )

    N = int(unique_N[0])
    L = N + 1

    # --------------------------------------------------------
    # Verify canonical start and target.
    # --------------------------------------------------------
    canonical_checks = {
        "start_x": 0.0,
        "start_y": 0.0,
        "start_z": 0.0,
        "target_x": 1.0,
        "target_y": 0.0,
        "target_z": 0.0,
    }

    for col, expected in canonical_checks.items():
        max_error = float(
            np.max(
                np.abs(
                    df[col].to_numpy(dtype=float)
                    - expected
                )
            )
        )

        if max_error > 1e-9:
            raise ValueError(
                "{} is not canonical. Maximum error = {}".format(
                    col,
                    max_error,
                )
            )

    # --------------------------------------------------------
    # Inputs
    # --------------------------------------------------------
    sphere_x = df["sphere_x"].to_numpy(dtype=np.float64)
    sphere_y = df["sphere_y"].to_numpy(dtype=np.float64)
    sphere_z = df["sphere_z"].to_numpy(dtype=np.float64)
    radius = df["sphere_radius"].to_numpy(dtype=np.float64)

    rho = np.sqrt(
        sphere_y * sphere_y
        + sphere_z * sphere_z
    )

    theta = np.arctan2(
        sphere_z,
        sphere_y,
    )

    X = np.column_stack(
        [
            sphere_x,
            rho,
            radius,
        ]
    )

    # --------------------------------------------------------
    # Outputs: full requested 4*(N+1) Bernstein coefficients
    #
    # [Px, Py_aligned, Pz_aligned, Ppsi]
    # --------------------------------------------------------
    px_cols = [
        "px_cp_{}".format(i)
        for i in range(L)
    ]

    py_cols = [
        "py_cp_{}".format(i)
        for i in range(L)
    ]

    pz_cols = [
        "pz_cp_{}".format(i)
        for i in range(L)
    ]

    psi_cols = [
        "psi_cp_{}".format(i)
        for i in range(L)
    ]

    for col in (
        px_cols
        + py_cols
        + pz_cols
        + psi_cols
    ):
        if col not in df.columns:
            raise ValueError(
                "Dataset is missing output column: {}".format(col)
            )

    px = df[px_cols].to_numpy(dtype=np.float64)
    py = df[py_cols].to_numpy(dtype=np.float64)
    pz = df[pz_cols].to_numpy(dtype=np.float64)
    psi = df[psi_cols].to_numpy(dtype=np.float64)

    theta_matrix = theta[:, None]

    py_aligned, pz_aligned = rotate_yz_to_obstacle_aligned(
        py,
        pz,
        theta_matrix,
    )

    Y = np.concatenate(
        [
            px,
            py_aligned,
            pz_aligned,
            psi,
        ],
        axis=1,
    )

    metadata = {
        "N": N,
        "L": L,
        "theta": theta,
        "rho": rho,
        "case_type": df["case_type"].astype(str).to_numpy(),
        "sample_id": df["sample_id"].to_numpy(),
        "sphere_x": sphere_x,
        "sphere_y": sphere_y,
        "sphere_z": sphere_z,
        "sphere_radius": radius,
    }

    return df, X, Y, metadata


# ============================================================
# Stratified train/validation/test split
#
# Keeps blocking / near_miss / clear proportions represented
# in all three subsets without sklearn.
# ============================================================
def stratified_split(
    case_types,
    train_fraction,
    val_fraction,
    seed,
):
    rng = np.random.default_rng(seed)

    train_indices = []
    val_indices = []
    test_indices = []

    unique_types = sorted(set(case_types.tolist()))

    for case_type in unique_types:
        indices = np.where(case_types == case_type)[0]

        rng.shuffle(indices)

        n = len(indices)

        n_train = int(round(train_fraction * n))
        n_val = int(round(val_fraction * n))

        # Guarantee at least one test point when possible.
        if n >= 3:
            n_train = min(n_train, n - 2)
            n_val = min(n_val, n - n_train - 1)

        train_indices.extend(
            indices[:n_train].tolist()
        )

        val_indices.extend(
            indices[
                n_train:
                n_train + n_val
            ].tolist()
        )

        test_indices.extend(
            indices[
                n_train + n_val:
            ].tolist()
        )

    rng.shuffle(train_indices)
    rng.shuffle(val_indices)
    rng.shuffle(test_indices)

    return (
        np.array(train_indices, dtype=np.int64),
        np.array(val_indices, dtype=np.int64),
        np.array(test_indices, dtype=np.int64),
    )


# ============================================================
# Standardization
# ============================================================
def compute_normalizer(array):
    mean = array.mean(axis=0)
    std = array.std(axis=0)

    # Constant or nearly constant columns (e.g. yaw in this dataset)
    # should not cause division by zero.
    std_safe = std.copy()
    std_safe[std_safe < 1e-10] = 1.0

    return mean, std_safe


def normalize(array, mean, std):
    return (array - mean) / std


def denormalize(array, mean, std):
    return array * std + mean


# ============================================================
# Neural network
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
# Training helpers
# ============================================================
def mse_on_loader(model, loader, criterion, device):
    model.eval()

    total_loss = 0.0
    total_count = 0

    with torch.no_grad():
        for xb, yb in loader:
            xb = xb.to(device)
            yb = yb.to(device)

            pred = model(xb)
            loss = criterion(pred, yb)

            batch_size = xb.shape[0]

            total_loss += (
                float(loss.item())
                * batch_size
            )

            total_count += batch_size

    return total_loss / max(total_count, 1)


def predict_array(model, X_norm, device, batch_size=1024):
    dataset = TensorDataset(
        torch.tensor(
            X_norm,
            dtype=torch.float32,
        )
    )

    loader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=False,
    )

    predictions = []

    model.eval()

    with torch.no_grad():
        for (xb,) in loader:
            xb = xb.to(device)

            pred = model(xb)

            predictions.append(
                pred.cpu().numpy()
            )

    return np.concatenate(
        predictions,
        axis=0,
    )


# ============================================================
# Dense Bernstein geometry evaluation
# ============================================================
def bernstein_basis_matrix(N, samples):
    s = np.linspace(
        0.0,
        1.0,
        samples,
    )

    basis = np.zeros(
        (samples, N + 1),
        dtype=np.float64,
    )

    for i in range(N + 1):
        basis[:, i] = (
            math.comb(N, i)
            * (s ** i)
            * ((1.0 - s) ** (N - i))
        )

    return basis


def unpack_outputs(Y, L):
    px = Y[:, 0:L]
    py = Y[:, L:2 * L]
    pz = Y[:, 2 * L:3 * L]
    psi = Y[:, 3 * L:4 * L]

    return px, py, pz, psi


def dense_clearances_aligned(
    Y_aligned,
    X,
    N,
    samples=250,
):
    L = N + 1

    px, py, pz, _ = unpack_outputs(
        Y_aligned,
        L,
    )

    B = bernstein_basis_matrix(
        N,
        samples,
    )

    # row trajectory:
    # [num_cases, samples]
    curve_x = px @ B.T
    curve_y = py @ B.T
    curve_z = pz @ B.T

    sphere_x = X[:, 0][:, None]
    sphere_rho = X[:, 1][:, None]
    radius = X[:, 2]

    distance = np.sqrt(
        (curve_x - sphere_x) ** 2
        + (curve_y - sphere_rho) ** 2
        + curve_z ** 2
    )

    minimum_distance = distance.min(axis=1)

    clearance = (
        minimum_distance
        - radius
    )

    return clearance


# ============================================================
# Diagnostics
# ============================================================
def group_metrics(
    Y_true,
    Y_pred,
    L,
):
    names = [
        ("Px", 0, L),
        ("Py_aligned", L, 2 * L),
        ("Pz_aligned", 2 * L, 3 * L),
        ("Psi", 3 * L, 4 * L),
    ]

    rows = []

    for name, a, b in names:
        err = (
            Y_pred[:, a:b]
            - Y_true[:, a:b]
        )

        rmse = float(
            np.sqrt(
                np.mean(err * err)
            )
        )

        mae = float(
            np.mean(
                np.abs(err)
            )
        )

        rows.append(
            (name, rmse, mae)
        )

    return rows


def save_learning_curve(
    train_history,
    val_history,
    output_path,
):
    fig, ax = plt.subplots(
        figsize=(9, 6)
    )

    epochs = np.arange(
        1,
        len(train_history) + 1,
    )

    ax.plot(
        epochs,
        train_history,
        label="Train loss",
    )

    ax.plot(
        epochs,
        val_history,
        label="Validation loss",
    )

    ax.set_xlabel("Epoch")
    ax.set_ylabel("Normalized MSE")
    ax.set_title(
        "Training history"
    )

    ax.set_yscale("log")
    ax.grid(alpha=0.25)
    ax.legend()

    fig.tight_layout()

    fig.savefig(
        output_path,
        dpi=220,
        bbox_inches="tight",
    )

    plt.close(fig)


def save_prediction_scatter(
    Y_true,
    Y_pred,
    output_path,
):
    true_flat = Y_true.reshape(-1)
    pred_flat = Y_pred.reshape(-1)

    fig, ax = plt.subplots(
        figsize=(7, 7)
    )

    ax.scatter(
        true_flat,
        pred_flat,
        s=10,
        alpha=0.35,
    )

    low = float(
        min(
            true_flat.min(),
            pred_flat.min(),
        )
    )

    high = float(
        max(
            true_flat.max(),
            pred_flat.max(),
        )
    )

    ax.plot(
        [low, high],
        [low, high],
        linestyle="--",
        linewidth=1.5,
        label="Perfect prediction",
    )

    ax.set_xlabel(
        "OCP Bernstein coefficient"
    )

    ax.set_ylabel(
        "NN-predicted coefficient"
    )

    ax.set_title(
        "Test-set coefficient predictions"
    )

    ax.grid(alpha=0.25)
    ax.legend()

    fig.tight_layout()

    fig.savefig(
        output_path,
        dpi=220,
        bbox_inches="tight",
    )

    plt.close(fig)


def save_clearance_plot(
    true_clearance,
    pred_clearance,
    case_types,
    output_path,
):
    fig, ax = plt.subplots(
        figsize=(9, 7)
    )

    for case_type in [
        "blocking",
        "near_miss",
        "clear",
    ]:
        mask = (
            case_types
            == case_type
        )

        if not np.any(mask):
            continue

        ax.scatter(
            true_clearance[mask],
            pred_clearance[mask],
            s=28,
            alpha=0.6,
            label=case_type,
        )

    low = float(
        min(
            true_clearance.min(),
            pred_clearance.min(),
        )
    )

    high = float(
        max(
            true_clearance.max(),
            pred_clearance.max(),
        )
    )

    ax.plot(
        [low, high],
        [low, high],
        linestyle="--",
        linewidth=1.5,
        label="Perfect agreement",
    )

    ax.axhline(
        0.0,
        linestyle=":",
        linewidth=1.5,
        label="Predicted collision boundary",
    )

    ax.set_xlabel(
        "OCP trajectory clearance"
    )

    ax.set_ylabel(
        "NN trajectory clearance"
    )

    ax.set_title(
        "Test-set sphere clearance"
    )

    ax.grid(alpha=0.25)
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
            "Train a CPU-only PyTorch residual MLP to predict "
            "Bernstein control points for the canonical "
            "single-spherical-obstacle OCP."
        )
    )

    parser.add_argument(
        "dataset",
        nargs="?",
        default=(
            "single_drone_single_obs_training_dataset_001.csv"
        ),
        help="Training CSV",
    )

    parser.add_argument(
        "--epochs",
        type=int,
        default=500,
    )

    parser.add_argument(
        "--batch-size",
        type=int,
        default=256,
    )

    parser.add_argument(
        "--lr",
        type=float,
        default=1e-3,
    )

    parser.add_argument(
        "--weight-decay",
        type=float,
        default=1e-5,
    )

    parser.add_argument(
        "--patience",
        type=int,
        default=60,
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=20260821,
    )

    parser.add_argument(
        "--output-prefix",
        default="single_obs_mlp",
    )

    args = parser.parse_args()

    set_seed(args.seed)

    # CPU only, intentionally.
    device = torch.device("cpu")

    dataset_path = Path(args.dataset)

    if not dataset_path.exists():
        raise FileNotFoundError(
            "Dataset not found: {}".format(
                dataset_path
            )
        )

    print(
        "============================================================"
    )
    print(
        "CANONICAL SINGLE-SPHERE BERNSTEIN NN TRAINING"
    )
    print(
        "============================================================"
    )
    print(
        "PyTorch version : {}".format(
            torch.__version__
        )
    )
    print(
        "Device          : {}".format(
            device
        )
    )
    print(
        "Dataset         : {}".format(
            dataset_path
        )
    )

    df, X, Y, metadata = prepare_dataset(
        dataset_path
    )

    N = metadata["N"]
    L = metadata["L"]

    print(
        "Samples         : {}".format(
            len(df)
        )
    )
    print(
        "Bernstein N     : {}".format(
            N
        )
    )
    print(
        "NN inputs       : 3 "
        "[sphere_x, rho, sphere_radius]"
    )
    print(
        "NN outputs      : {} "
        "[Px, Py_aligned, Pz_aligned, Psi]".format(
            4 * L
        )
    )

    print("\nCase distribution:")
    print(
        df["case_type"]
        .value_counts()
        .to_string()
    )

    # Useful diagnostic before training.
    psi_part = Y[:, 3 * L:4 * L]

    print(
        "\nMaximum absolute Psi control point "
        "in dataset: {:.12e}".format(
            float(
                np.max(
                    np.abs(psi_part)
                )
            )
        )
    )

    # --------------------------------------------------------
    # Split
    # --------------------------------------------------------
    train_idx, val_idx, test_idx = (
        stratified_split(
            metadata["case_type"],
            train_fraction=0.80,
            val_fraction=0.10,
            seed=args.seed,
        )
    )

    print(
        "\nSplit:"
    )
    print(
        "  train = {}".format(
            len(train_idx)
        )
    )
    print(
        "  val   = {}".format(
            len(val_idx)
        )
    )
    print(
        "  test  = {}".format(
            len(test_idx)
        )
    )

    X_train = X[train_idx]
    X_val = X[val_idx]
    X_test = X[test_idx]

    Y_train = Y[train_idx]
    Y_val = Y[val_idx]
    Y_test = Y[test_idx]

    # --------------------------------------------------------
    # Normalize using TRAIN ONLY.
    # --------------------------------------------------------
    X_mean, X_std = compute_normalizer(
        X_train
    )

    Y_mean, Y_std = compute_normalizer(
        Y_train
    )

    X_train_n = normalize(
        X_train,
        X_mean,
        X_std,
    )

    X_val_n = normalize(
        X_val,
        X_mean,
        X_std,
    )

    X_test_n = normalize(
        X_test,
        X_mean,
        X_std,
    )

    Y_train_n = normalize(
        Y_train,
        Y_mean,
        Y_std,
    )

    Y_val_n = normalize(
        Y_val,
        Y_mean,
        Y_std,
    )

    # --------------------------------------------------------
    # DataLoaders
    # --------------------------------------------------------
    train_dataset = TensorDataset(
        torch.tensor(
            X_train_n,
            dtype=torch.float32,
        ),
        torch.tensor(
            Y_train_n,
            dtype=torch.float32,
        ),
    )

    val_dataset = TensorDataset(
        torch.tensor(
            X_val_n,
            dtype=torch.float32,
        ),
        torch.tensor(
            Y_val_n,
            dtype=torch.float32,
        ),
    )

    train_loader = DataLoader(
        train_dataset,
        batch_size=args.batch_size,
        shuffle=True,
        num_workers=0,
    )

    val_loader = DataLoader(
        val_dataset,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=0,
    )

    # --------------------------------------------------------
    # Model
    # --------------------------------------------------------
    model = BernsteinMLP(
        input_dim=X.shape[1],
        output_dim=Y.shape[1],
        width=128,
        residual_blocks=3,
    ).to(device)

    number_of_parameters = sum(
        p.numel()
        for p in model.parameters()
        if p.requires_grad
    )

    print(
        "\nTrainable parameters: {}".format(
            number_of_parameters
        )
    )

    criterion = nn.MSELoss()

    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=args.lr,
        weight_decay=args.weight_decay,
    )

    # Reduce learning rate when validation progress stalls.
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode="min",
        factor=0.5,
        patience=15,
        threshold=1e-5,
        min_lr=1e-6,
    )

    # --------------------------------------------------------
    # Training
    # --------------------------------------------------------
    train_history = []
    val_history = []

    best_val = float("inf")
    best_state = None
    epochs_without_improvement = 0
    best_epoch = 0

    print(
        "\nStarting training...\n"
    )

    for epoch in range(1, args.epochs + 1):
        model.train()

        running_loss = 0.0
        running_count = 0

        for xb, yb in train_loader:
            xb = xb.to(device)
            yb = yb.to(device)

            optimizer.zero_grad(
                set_to_none=True
            )

            prediction = model(xb)

            loss = criterion(
                prediction,
                yb,
            )

            loss.backward()

            optimizer.step()

            batch_size = xb.shape[0]

            running_loss += (
                float(loss.item())
                * batch_size
            )

            running_count += batch_size

        train_loss = (
            running_loss
            / max(running_count, 1)
        )

        val_loss = mse_on_loader(
            model,
            val_loader,
            criterion,
            device,
        )

        train_history.append(
            train_loss
        )

        val_history.append(
            val_loss
        )

        scheduler.step(
            val_loss
        )

        improved = (
            val_loss
            < best_val - 1e-8
        )

        if improved:
            best_val = val_loss
            best_epoch = epoch

            best_state = {
                key: value.detach().cpu().clone()
                for key, value
                in model.state_dict().items()
            }

            epochs_without_improvement = 0
        else:
            epochs_without_improvement += 1

        if (
            epoch == 1
            or epoch % 10 == 0
            or improved
        ):
            current_lr = (
                optimizer.param_groups[0]["lr"]
            )

            print(
                "Epoch {:4d} | "
                "train {:.8e} | "
                "val {:.8e} | "
                "lr {:.3e}{}".format(
                    epoch,
                    train_loss,
                    val_loss,
                    current_lr,
                    "  *" if improved else "",
                )
            )

        if (
            epochs_without_improvement
            >= args.patience
        ):
            print(
                "\nEarly stopping at epoch {}. "
                "Best validation epoch = {}.".format(
                    epoch,
                    best_epoch,
                )
            )

            break

    if best_state is None:
        raise RuntimeError(
            "Training did not produce a valid model state."
        )

    model.load_state_dict(
        best_state
    )

    # --------------------------------------------------------
    # Test prediction in original canonical coefficient units
    # --------------------------------------------------------
    pred_test_n = predict_array(
        model,
        X_test_n,
        device,
    )

    pred_test = denormalize(
        pred_test_n,
        Y_mean,
        Y_std,
    )

    error = (
        pred_test
        - Y_test
    )

    overall_rmse = float(
        np.sqrt(
            np.mean(
                error * error
            )
        )
    )

    overall_mae = float(
        np.mean(
            np.abs(error)
        )
    )

    print(
        "\n============================================================"
    )
    print(
        "TEST RESULTS"
    )
    print(
        "============================================================"
    )

    print(
        "Overall coefficient RMSE : {:.10e}".format(
            overall_rmse
        )
    )

    print(
        "Overall coefficient MAE  : {:.10e}".format(
            overall_mae
        )
    )

    print(
        "\nPer-output-group errors:"
    )

    metrics = group_metrics(
        Y_test,
        pred_test,
        L,
    )

    for name, rmse, mae in metrics:
        print(
            "  {:12s} RMSE = {:.10e} | MAE = {:.10e}".format(
                name,
                rmse,
                mae,
            )
        )

    # --------------------------------------------------------
    # Geometric safety diagnostic
    #
    # This is NOT part of the training loss yet.
    # It tells us how often pure supervised regression produces
    # a dense reconstructed path that intersects the sphere.
    # --------------------------------------------------------
    true_clearance = dense_clearances_aligned(
        Y_test,
        X_test,
        N,
        samples=300,
    )

    pred_clearance = dense_clearances_aligned(
        pred_test,
        X_test,
        N,
        samples=300,
    )

    predicted_collisions_0 = int(
        np.sum(
            pred_clearance < 0.0
        )
    )

    predicted_collisions_5mm = int(
        np.sum(
            pred_clearance < -0.005
        )
    )

    print(
        "\nDense reconstructed trajectory check:"
    )

    print(
        "  Test OCP min clearance       : {:.10e}".format(
            float(
                true_clearance.min()
            )
        )
    )

    print(
        "  NN predicted min clearance   : {:.10e}".format(
            float(
                pred_clearance.min()
            )
        )
    )

    print(
        "  NN cases with clearance < 0  : {} / {}".format(
            predicted_collisions_0,
            len(test_idx),
        )
    )

    print(
        "  NN penetration > 0.005       : {} / {}".format(
            predicted_collisions_5mm,
            len(test_idx),
        )
    )

    # --------------------------------------------------------
    # Save model + preprocessing
    # --------------------------------------------------------
    model_path = Path(
        args.output_prefix
        + "_best.pt"
    )

    checkpoint = {
        "model_state_dict": (
            model.state_dict()
        ),
        "input_dim": int(X.shape[1]),
        "output_dim": int(Y.shape[1]),
        "width": 128,
        "residual_blocks": 3,
        "N": int(N),
        "L": int(L),
        "X_mean": X_mean.astype(
            np.float64
        ),
        "X_std": X_std.astype(
            np.float64
        ),
        "Y_mean": Y_mean.astype(
            np.float64
        ),
        "Y_std": Y_std.astype(
            np.float64
        ),
        "input_definition": [
            "sphere_x",
            "sphere_rho",
            "sphere_radius",
        ],
        "output_definition": (
            ["px_cp_{}".format(i) for i in range(L)]
            + ["py_aligned_cp_{}".format(i) for i in range(L)]
            + ["pz_aligned_cp_{}".format(i) for i in range(L)]
            + ["psi_cp_{}".format(i) for i in range(L)]
        ),
        "geometry": (
            "Canonical start=[0,0,0], target=[1,0,0]. "
            "Obstacle and output y/z are rotated about x "
            "so obstacle lies at [sphere_x, rho, 0]."
        ),
        "seed": int(args.seed),
        "best_epoch": int(best_epoch),
        "best_validation_loss": float(
            best_val
        ),
        "test_coefficient_rmse": float(
            overall_rmse
        ),
        "test_coefficient_mae": float(
            overall_mae
        ),
        "test_predicted_collisions": int(
            predicted_collisions_0
        ),
    }

    torch.save(
        checkpoint,
        model_path,
    )

    # Human-readable normalization/settings file.
    metadata_path = Path(
        args.output_prefix
        + "_metadata.json"
    )

    json_data = {
        "dataset": str(dataset_path),
        "N": int(N),
        "L": int(L),
        "input_dim": int(X.shape[1]),
        "output_dim": int(Y.shape[1]),
        "input_definition": checkpoint[
            "input_definition"
        ],
        "output_definition": checkpoint[
            "output_definition"
        ],
        "X_mean": X_mean.tolist(),
        "X_std": X_std.tolist(),
        "Y_mean": Y_mean.tolist(),
        "Y_std": Y_std.tolist(),
        "best_epoch": int(best_epoch),
        "best_validation_loss": float(
            best_val
        ),
        "test_coefficient_rmse": float(
            overall_rmse
        ),
        "test_coefficient_mae": float(
            overall_mae
        ),
        "predicted_collision_count": int(
            predicted_collisions_0
        ),
        "test_count": int(
            len(test_idx)
        ),
    }

    with open(
        metadata_path,
        "w",
    ) as f:
        json.dump(
            json_data,
            f,
            indent=2,
        )

    # --------------------------------------------------------
    # Save plots
    # --------------------------------------------------------
    learning_curve_path = Path(
        args.output_prefix
        + "_learning_curve.png"
    )

    prediction_scatter_path = Path(
        args.output_prefix
        + "_prediction_scatter.png"
    )

    clearance_plot_path = Path(
        args.output_prefix
        + "_clearance_test.png"
    )

    save_learning_curve(
        train_history,
        val_history,
        learning_curve_path,
    )

    save_prediction_scatter(
        Y_test,
        pred_test,
        prediction_scatter_path,
    )

    save_clearance_plot(
        true_clearance,
        pred_clearance,
        metadata["case_type"][test_idx],
        clearance_plot_path,
    )

    print(
        "\nSaved model      : {}".format(
            model_path
        )
    )
    print(
        "Saved metadata   : {}".format(
            metadata_path
        )
    )
    print(
        "Saved curve plot : {}".format(
            learning_curve_path
        )
    )
    print(
        "Saved pred plot  : {}".format(
            prediction_scatter_path
        )
    )
    print(
        "Saved safety plot: {}".format(
            clearance_plot_path
        )
    )

    print(
        "============================================================"
    )


if __name__ == "__main__":
    main()
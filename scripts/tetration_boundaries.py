from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


E_NEG_E = float(np.exp(-np.e))
E_POS_INV_E = float(np.exp(1.0 / np.e))
E_NEG_INV_E = float(np.exp(-1.0 / np.e))


def c_cos_collatz(z: complex) -> complex:
    """Cosine-interpolated Collatz gate."""
    cz = np.cos(np.pi * z)
    return (z / 4.0) * (1.0 + cz) + ((3.0 * z + 1.0) / 2.0) * (1.0 - cz)


def blended_gate(z: complex, collatz_weight: float) -> complex:
    """
    Controlled switch from pure tetration to full Collatz gating:
      collatz_weight = 1 -> C(z)
      collatz_weight = 0 -> z
    """
    return collatz_weight * c_cos_collatz(z) + (1.0 - collatz_weight) * z


@dataclass
class OrbitResult:
    status: str  # converged | bounded | escaped
    final_z: complex
    max_abs: float
    tail_max_delta: float


def iterate_orbit(
    c: complex,
    *,
    collatz_weight: float,
    z0: complex,
    steps: int,
    radius: float,
    tol: float,
    tail_window: int,
    imag_tol: float,
    alpha: float,
    tanh_scale: float,
) -> OrbitResult:
    z = complex(z0)
    max_abs = abs(z)
    deltas: list[float] = []
    log_radius = float(np.log(radius))

    for _ in range(steps):
        gz = blended_gate(z, collatz_weight=collatz_weight)

        # Optional phase kick; keep alpha=0 for pure tetration-style scans.
        phi = float(alpha * np.tanh(np.real(gz) / tanh_scale))
        c_eff = c * np.exp(1j * phi)
        arg = gz * np.log(c_eff)

        if (not np.isfinite(arg.real)) or (not np.isfinite(arg.imag)):
            return OrbitResult("escaped", z, max_abs, math.inf)
        if arg.real > log_radius:
            return OrbitResult("escaped", z, max_abs, math.inf)

        z_next = np.exp(arg)
        if (not np.isfinite(z_next.real)) or (not np.isfinite(z_next.imag)):
            return OrbitResult("escaped", z, max_abs, math.inf)

        abs_next = abs(z_next)
        max_abs = max(max_abs, abs_next)
        if abs_next > radius:
            return OrbitResult("escaped", complex(z_next), max_abs, math.inf)

        delta = abs(z_next - z)
        deltas.append(delta)
        if len(deltas) > tail_window:
            deltas.pop(0)

        z = complex(z_next)

    tail_max_delta = max(deltas) if deltas else math.inf
    if len(deltas) == tail_window and tail_max_delta < tol and abs(np.imag(z)) < imag_tol:
        return OrbitResult("converged", z, max_abs, tail_max_delta)
    return OrbitResult("bounded", z, max_abs, tail_max_delta)


def parse_weights(raw: str) -> list[float]:
    weights = [float(x.strip()) for x in raw.split(",") if x.strip()]
    if not weights:
        raise ValueError("No collatz weights provided.")
    for w in weights:
        if (w < 0.0) or (w > 1.0):
            raise ValueError(f"collatz weight must be in [0,1], got {w}")
    return weights


def mask_to_intervals(xs: np.ndarray, mask: np.ndarray) -> list[tuple[float, float, int]]:
    intervals: list[tuple[float, float, int]] = []
    n = len(xs)
    i = 0
    while i < n:
        if not mask[i]:
            i += 1
            continue
        j = i
        while (j + 1) < n and mask[j + 1]:
            j += 1
        intervals.append((float(xs[i]), float(xs[j]), j - i + 1))
        i = j + 1
    return intervals


def select_primary_interval(intervals: list[tuple[float, float, int]]) -> tuple[float, float, int] | None:
    if not intervals:
        return None
    # Prefer most points; break ties with larger width.
    return max(intervals, key=lambda it: (it[2], it[1] - it[0]))


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Scan real-base tetration boundaries while smoothly turning off the "
            "Collatz gate via collatz_weight in [0,1]."
        )
    )
    parser.add_argument("--weights", default="1.0,0.8,0.6,0.4,0.2,0.0")
    parser.add_argument("--c-min", type=float, default=0.02)
    parser.add_argument("--c-max", type=float, default=1.90)
    parser.add_argument("--c-points", type=int, default=900)
    parser.add_argument("--z0", type=float, default=1.0, help="Initial real z0 for the orbit.")
    parser.add_argument("--steps", type=int, default=260)
    parser.add_argument("--radius", type=float, default=1e6)
    parser.add_argument("--tol", type=float, default=1e-9)
    parser.add_argument("--tail-window", type=int, default=20)
    parser.add_argument("--imag-tol", type=float, default=1e-8)
    parser.add_argument("--alpha", type=float, default=0.0, help="Optional phase-kick strength.")
    parser.add_argument("--tanh-scale", type=float, default=1.0)
    parser.add_argument("--scan-csv", default="outputs/csv/tetration_boundaries_scan.csv")
    parser.add_argument("--summary-csv", default="outputs/csv/tetration_boundaries_summary.csv")
    parser.add_argument("--plot-out", default="outputs/figures/tetration_boundaries_edges.png")
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()

    weights = parse_weights(args.weights)
    c_values = np.linspace(args.c_min, args.c_max, args.c_points, dtype=float)

    scan_csv = Path(args.scan_csv)
    summary_csv = Path(args.summary_csv)
    plot_out = Path(args.plot_out)
    ensure_parent(scan_csv)
    ensure_parent(summary_csv)
    if not args.no_plot:
        ensure_parent(plot_out)

    summary_rows: list[dict[str, float | int]] = []

    with scan_csv.open("w", newline="", encoding="utf-8") as f_scan:
        w_scan = csv.writer(f_scan)
        w_scan.writerow(
            [
                "collatz_weight",
                "c_real",
                "status",
                "final_real",
                "final_imag",
                "max_abs",
                "tail_max_delta",
            ]
        )

        for weight in weights:
            statuses: list[str] = []
            for c in c_values:
                res = iterate_orbit(
                    complex(c, 0.0),
                    collatz_weight=weight,
                    z0=complex(args.z0, 0.0),
                    steps=args.steps,
                    radius=args.radius,
                    tol=args.tol,
                    tail_window=args.tail_window,
                    imag_tol=args.imag_tol,
                    alpha=args.alpha,
                    tanh_scale=args.tanh_scale,
                )
                statuses.append(res.status)
                w_scan.writerow(
                    [
                        f"{weight:.6f}",
                        f"{c:.12g}",
                        res.status,
                        f"{res.final_z.real:.12g}",
                        f"{res.final_z.imag:.12g}",
                        f"{res.max_abs:.12g}",
                        f"{res.tail_max_delta:.12g}",
                    ]
                )

            conv_mask = np.array([s == "converged" for s in statuses], dtype=bool)
            intervals = mask_to_intervals(c_values, conv_mask)
            primary = select_primary_interval(intervals)

            if primary is None:
                lo = float("nan")
                hi = float("nan")
                width = 0.0
                n_points = 0
            else:
                lo, hi, n_points = primary
                width = hi - lo

            summary_rows.append(
                {
                    "collatz_weight": weight,
                    "interval_lo": lo,
                    "interval_hi": hi,
                    "interval_width": width,
                    "interval_points": int(n_points),
                    "extend_below_exp_neg_e": max(0.0, E_NEG_E - lo) if np.isfinite(lo) else float("nan"),
                    "extend_above_exp_pos_inv_e": max(0.0, hi - E_POS_INV_E) if np.isfinite(hi) else float("nan"),
                    "extend_above_exp_neg_inv_e": max(0.0, hi - E_NEG_INV_E) if np.isfinite(hi) else float("nan"),
                }
            )

    with summary_csv.open("w", newline="", encoding="utf-8") as f_sum:
        w_sum = csv.writer(f_sum)
        w_sum.writerow(
            [
                "collatz_weight",
                "interval_lo",
                "interval_hi",
                "interval_width",
                "interval_points",
                "extend_below_exp_neg_e",
                "extend_above_exp_pos_inv_e",
                "extend_above_exp_neg_inv_e",
            ]
        )
        for row in summary_rows:
            w_sum.writerow(
                [
                    f"{row['collatz_weight']:.6f}",
                    f"{row['interval_lo']:.12g}" if np.isfinite(row["interval_lo"]) else "nan",
                    f"{row['interval_hi']:.12g}" if np.isfinite(row["interval_hi"]) else "nan",
                    f"{row['interval_width']:.12g}",
                    f"{int(row['interval_points'])}",
                    f"{row['extend_below_exp_neg_e']:.12g}" if np.isfinite(row["extend_below_exp_neg_e"]) else "nan",
                    f"{row['extend_above_exp_pos_inv_e']:.12g}" if np.isfinite(row["extend_above_exp_pos_inv_e"]) else "nan",
                    f"{row['extend_above_exp_neg_inv_e']:.12g}" if np.isfinite(row["extend_above_exp_neg_inv_e"]) else "nan",
                ]
            )

    print("Reference boundaries:")
    print(f"  exp(-e)   = {E_NEG_E:.12g}")
    print(f"  exp(1/e)  = {E_POS_INV_E:.12g}  (classical upper bound)")
    print(f"  exp(-1/e) = {E_NEG_INV_E:.12g}  (included for comparison)")
    print("")
    print("Primary convergent interval per collatz_weight:")
    print("  weight    lo            hi            width         ext_below_e^-e  ext_above_e^(1/e)")
    for row in summary_rows:
        lo = row["interval_lo"]
        hi = row["interval_hi"]
        width = row["interval_width"]
        ext_lo = row["extend_below_exp_neg_e"]
        ext_hi = row["extend_above_exp_pos_inv_e"]
        print(
            f"  {row['collatz_weight']:>6.3f}   "
            f"{lo:>11.6f}   {hi:>11.6f}   {width:>11.6f}   "
            f"{ext_lo:>13.6f}   {ext_hi:>13.6f}"
        )

    if not args.no_plot:
        xs = np.array([float(row["collatz_weight"]) for row in summary_rows], dtype=float)
        ys_lo = np.array([float(row["interval_lo"]) for row in summary_rows], dtype=float)
        ys_hi = np.array([float(row["interval_hi"]) for row in summary_rows], dtype=float)
        ys_width = np.array([float(row["interval_width"]) for row in summary_rows], dtype=float)

        fig, axes = plt.subplots(2, 1, figsize=(8.6, 8.4), sharex=True)
        ax0, ax1 = axes

        ax0.plot(xs, ys_lo, marker="o", label="primary interval lo")
        ax0.plot(xs, ys_hi, marker="o", label="primary interval hi")
        ax0.axhline(E_NEG_E, linestyle="--", color="tab:green", label="exp(-e)")
        ax0.axhline(E_POS_INV_E, linestyle="--", color="tab:red", label="exp(1/e)")
        ax0.axhline(E_NEG_INV_E, linestyle=":", color="tab:purple", label="exp(-1/e)")
        ax0.set_ylabel("c (real)")
        ax0.set_title("Convergent real-base interval vs collatz_weight")
        ax0.grid(alpha=0.25)
        ax0.legend(loc="best")

        ax1.plot(xs, ys_width, marker="o", color="tab:blue")
        ax1.set_xlabel("collatz_weight (1 = full Collatz, 0 = pure tetration)")
        ax1.set_ylabel("primary interval width")
        ax1.grid(alpha=0.25)

        fig.tight_layout()
        fig.savefig(plot_out, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"\nSaved plot: {plot_out}")

    print(f"Saved scan CSV: {scan_csv}")
    print(f"Saved summary CSV: {summary_csv}")


if __name__ == "__main__":
    main()

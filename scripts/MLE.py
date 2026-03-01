from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt

# Allow running directly from repo root without package install.
REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from collatzscape.metrics import finite_time_lyapunov_rough, maximal_lyapunov_stats


def _parse_float_list(csv_values: str) -> list[float]:
    return [float(v.strip()) for v in csv_values.split(",") if v.strip()]


def _parse_int_list(csv_values: str) -> list[int]:
    return [int(v.strip()) for v in csv_values.split(",") if v.strip()]


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Accurate MLE sweep for Collatzscape (rough vs Benettin-style)."
    )
    ap.add_argument("--c-re", type=float, default=1.8725)
    ap.add_argument("--c-im", type=float, default=0.0)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--steps", type=int, default=120)
    ap.add_argument("--gammas", type=str, default="0.0,0.05,0.1,0.15,0.3,0.5")
    ap.add_argument("--directions", type=str, default="1,-1")
    ap.add_argument("--kick-mode", type=str, default="tanh")
    ap.add_argument("--C-mode", type=str, default="cosine")
    ap.add_argument("--unpaired", action="store_true", help="Use w=conj(z) each step.")
    ap.add_argument("--z0-re", type=float, default=1.0)
    ap.add_argument("--z0-im", type=float, default=0.0)
    ap.add_argument("--eps", type=float, default=1e-8)
    ap.add_argument("--mle-trials", type=int, default=8)
    ap.add_argument("--mle-warmup-steps", type=int, default=32)
    ap.add_argument("--rough-warmup-steps", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-csv", type=str, default="outputs/mle_sweep.csv")
    ap.add_argument("--out-plot", type=str, default="figures/mle_sweep.png")
    args = ap.parse_args()

    c = complex(args.c_re, args.c_im)
    z0 = complex(args.z0_re, args.z0_im)
    paired = not args.unpaired
    gammas = _parse_float_list(args.gammas)
    directions = _parse_int_list(args.directions)

    out_csv = Path(args.out_csv)
    out_plot = Path(args.out_plot)
    os.makedirs(out_csv.parent, exist_ok=True)
    os.makedirs(out_plot.parent, exist_ok=True)

    rows: list[tuple[float, int, float, float, float, int]] = []
    for dir_idx, direction in enumerate(directions):
        for gamma_idx, gamma in enumerate(gammas):
            rough = finite_time_lyapunov_rough(
                z0,
                c=c,
                alpha=args.alpha,
                gamma=gamma,
                direction=direction,
                steps=args.steps,
                kick_mode=args.kick_mode,
                C_mode=args.C_mode,
                paired=paired,
                eps=args.eps,
                warmup_steps=args.rough_warmup_steps,
                seed=args.seed + 20000 + 1000 * dir_idx + gamma_idx,
            )
            mle_stats = maximal_lyapunov_stats(
                z0,
                c=c,
                alpha=args.alpha,
                gamma=gamma,
                direction=direction,
                steps=args.steps,
                kick_mode=args.kick_mode,
                C_mode=args.C_mode,
                paired=paired,
                eps=args.eps,
                warmup_steps=args.mle_warmup_steps,
                n_trials=args.mle_trials,
                seed=args.seed + 1000 * dir_idx + gamma_idx,
            )
            rows.append(
                (
                    gamma,
                    direction,
                    float(rough),
                    float(mle_stats["mle"]),
                    float(mle_stats["mle_std"]),
                    int(mle_stats["valid_trials"]),
                )
            )

    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["gamma", "direction", "ft_lyapunov_rough", "mle", "mle_std", "mle_valid_trials"])
        w.writerows(rows)

    fig, axes = plt.subplots(1, len(directions), figsize=(6.0 * max(1, len(directions)), 4.4), sharey=True)
    if len(directions) == 1:
        axes = [axes]

    for ax, direction in zip(axes, directions):
        dir_rows = sorted((r for r in rows if r[1] == direction), key=lambda r: r[0])
        xs = [r[0] for r in dir_rows]
        rough = [r[2] for r in dir_rows]
        mle = [r[3] for r in dir_rows]
        mle_std = [r[4] for r in dir_rows]
        ax.plot(xs, rough, marker="o", label="rough")
        ax.errorbar(xs, mle, yerr=mle_std, marker="o", capsize=3, label="accurate MLE")
        ax.set_title(f"direction={direction}")
        ax.set_xlabel("gamma")
        ax.grid(alpha=0.25)

    axes[0].set_ylabel("Lyapunov exponent (per step)")
    axes[0].legend()
    fig.suptitle("Collatzscape MLE: rough vs accurate estimator")
    fig.tight_layout()
    fig.savefig(out_plot, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print("Saved:")
    print(f"  csv : {out_csv.resolve()}")
    print(f"  plot: {out_plot.resolve()}")


if __name__ == "__main__":
    main()

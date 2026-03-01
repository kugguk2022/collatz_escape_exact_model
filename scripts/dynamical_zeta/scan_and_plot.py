from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from find_cycles import CycleRecord, discover_cycles, write_cycle_table
from finite_product import CycleData, evaluate_curve


def parse_complex(raw: str) -> complex:
    return complex(raw.strip().replace("i", "j"))


def parse_complex_list(raw: str) -> list[complex]:
    parts = [p.strip() for p in raw.split(",") if p.strip()]
    if not parts:
        raise ValueError("No c values provided")
    return [parse_complex(p) for p in parts]


def slug_complex(c: complex) -> str:
    r = f"{c.real:+.6f}".replace("+", "p").replace("-", "m")
    i = f"{c.imag:+.6f}".replace("+", "p").replace("-", "m")
    return f"c_{r}_{i}"


def cycles_to_cycledata(cycles: list[CycleRecord]) -> list[CycleData]:
    out: list[CycleData] = []
    for c in cycles:
        out.append(
            CycleData(
                period=c.period,
                multiplier=c.multiplier,
                location=c.location,
                hits=c.hits,
            )
        )
    return out


def write_combined_curve_csv(path: Path, rows: list[tuple[complex, float, complex]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["c_real", "c_imag", "s", "zeta_real", "zeta_imag", "abs_zeta", "log10_abs_zeta"])
        for c, s, z in rows:
            az = abs(z)
            w.writerow(
                [
                    f"{c.real:.16g}",
                    f"{c.imag:.16g}",
                    f"{s:.16g}",
                    f"{z.real:.16g}",
                    f"{z.imag:.16g}",
                    f"{az:.16g}",
                    f"{np.log10(max(az, 1e-300)):.16g}",
                ]
            )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "For several c values (including bay-boundary-like defaults), "
            "find cycles and plot |zeta(s)| from finite products."
        )
    )
    parser.add_argument(
        "--c-values",
        default="1.8725+0j,1.93+0j,1.96+0j",
        help="Comma-separated list of complex c values.",
    )
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--gamma", type=float, default=1.0)
    parser.add_argument("--direction", type=int, choices=[-1, 1], default=1)
    parser.add_argument("--kick-mode", choices=["tanh", "clip"], default="tanh")
    parser.add_argument("--C-mode", choices=["cosine", "linear"], default="cosine")

    parser.add_argument("--re-min", type=float, default=-2.0)
    parser.add_argument("--re-max", type=float, default=2.0)
    parser.add_argument("--im-min", type=float, default=-2.0)
    parser.add_argument("--im-max", type=float, default=2.0)
    parser.add_argument("--n-re", type=int, default=100)
    parser.add_argument("--n-im", type=int, default=100)
    parser.add_argument("--steps", type=int, default=220)
    parser.add_argument("--burn-in", type=int, default=80)
    parser.add_argument("--max-period", type=int, default=24)
    parser.add_argument("--confirm-depth", type=int, default=4)
    parser.add_argument("--close-tol", type=float, default=1e-5)
    parser.add_argument("--dedup-tol", type=float, default=3e-4)
    parser.add_argument("--escape-radius", type=float, default=1e6)
    parser.add_argument("--max-cycles", type=int, default=2000)

    parser.add_argument("--s-min", type=float, default=0.0)
    parser.add_argument("--s-max", type=float, default=3.0)
    parser.add_argument("--s-points", type=int, default=600)
    parser.add_argument("--multiplier-floor", type=float, default=1e-6)
    parser.add_argument("--product-max-period", type=int, default=0, help="0 means no cutoff.")
    parser.add_argument("--use-multiplier-abs", action="store_true")

    parser.add_argument("--cycles-dir", default="outputs/csv/dynamical_zeta")
    parser.add_argument("--curve-csv", default="outputs/csv/dynamical_zeta/multi_c_zeta_curves.csv")
    parser.add_argument("--plot-out", default="outputs/figures/dynamical_zeta_multi_c_abs_zeta.png")
    parser.add_argument("--log-y", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    c_values = parse_complex_list(args.c_values)
    cycles_dir = Path(args.cycles_dir).resolve()
    cycles_dir.mkdir(parents=True, exist_ok=True)
    max_period = None if int(args.product_max_period) <= 0 else int(args.product_max_period)

    fig = plt.figure(figsize=(9.0, 5.4))
    all_rows: list[tuple[complex, float, complex]] = []

    for c in c_values:
        cycles, stats = discover_cycles(
            c=c,
            alpha=args.alpha,
            gamma=args.gamma,
            direction=args.direction,
            kick_mode=args.kick_mode,
            C_mode=args.C_mode,
            re_min=args.re_min,
            re_max=args.re_max,
            im_min=args.im_min,
            im_max=args.im_max,
            n_re=args.n_re,
            n_im=args.n_im,
            steps=args.steps,
            burn_in=args.burn_in,
            max_period=args.max_period,
            confirm_depth=args.confirm_depth,
            close_tol=args.close_tol,
            dedup_tol=args.dedup_tol,
            escape_radius=args.escape_radius,
            max_cycles=args.max_cycles,
            verbose=args.verbose,
        )
        cycle_file = cycles_dir / f"cycles_{slug_complex(c)}.txt"
        write_cycle_table(
            cycle_file,
            cycles,
            c=c,
            alpha=args.alpha,
            gamma=args.gamma,
            direction=args.direction,
            kick_mode=args.kick_mode,
            C_mode=args.C_mode,
            stats=stats,
        )
        if not cycles:
            if args.verbose:
                print(f"[scan_and_plot] no cycles detected for c={c}, skipping zeta curve")
            continue

        cycle_data = cycles_to_cycledata(cycles)
        s_values, zeta = evaluate_curve(
            cycle_data,
            s_min=args.s_min,
            s_max=args.s_max,
            s_points=args.s_points,
            multiplier_floor=args.multiplier_floor,
            max_period=max_period,
            use_multiplier_abs=bool(args.use_multiplier_abs),
        )

        label = f"c={c.real:.6g}{c.imag:+.6g}i, cycles={len(cycles)}"
        plt.plot(s_values, np.abs(zeta), lw=1.7, label=label)
        for s, z in zip(s_values, zeta):
            all_rows.append((c, float(s), complex(z)))

    plt.xlabel("s")
    plt.ylabel("|zeta(s)|")
    plt.title("Finite-product dynamical zeta from detected cycles")
    if args.log_y:
        plt.yscale("log")
    plt.grid(alpha=0.25)
    if plt.gca().lines:
        plt.legend(fontsize=8)
    plt.tight_layout()

    plot_out = Path(args.plot_out).resolve()
    plot_out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_out, dpi=180, bbox_inches="tight")
    plt.close(fig)

    curve_csv = Path(args.curve_csv).resolve()
    write_combined_curve_csv(curve_csv, all_rows)

    print("Dynamical zeta multi-c scan")
    print(f"  c values      : {len(c_values)}")
    print(f"  cycle dir     : {cycles_dir}")
    print(f"  curves csv    : {curve_csv}")
    print(f"  plot          : {plot_out}")


if __name__ == "__main__":
    main()

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


@dataclass
class CycleData:
    period: int
    multiplier: complex
    location: complex
    hits: int


def read_cycle_table(path: Path) -> list[CycleData]:
    rows: list[CycleData] = []
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader((line for line in f if not line.startswith("#")), delimiter="\t")
        for row in reader:
            rows.append(
                CycleData(
                    period=int(row["period"]),
                    multiplier=complex(float(row["multiplier_real"]), float(row["multiplier_imag"])),
                    location=complex(float(row["location_real"]), float(row["location_imag"])),
                    hits=int(row.get("hits", 1)),
                )
            )
    return rows


def finite_product_zeta(
    s: complex,
    cycles: list[CycleData],
    *,
    multiplier_floor: float,
    max_period: int | None,
    use_multiplier_abs: bool,
) -> complex:
    prod = 1.0 + 0.0j
    for rec in cycles:
        if max_period is not None and rec.period > max_period:
            continue
        lam = rec.multiplier
        lam_abs = float(abs(lam)) if np.isfinite(lam.real) and np.isfinite(lam.imag) else float("nan")
        if use_multiplier_abs:
            if not np.isfinite(lam_abs):
                lam_eff = float(multiplier_floor)
            else:
                lam_eff = max(lam_abs, multiplier_floor)
        else:
            if (not np.isfinite(lam.real)) or (not np.isfinite(lam.imag)) or (not np.isfinite(lam_abs)):
                lam_eff = complex(multiplier_floor, 0.0)
            else:
                lam_eff = lam if lam_abs >= multiplier_floor else complex(multiplier_floor, 0.0)
        term = np.exp(-s * rec.period) / lam_eff
        factor = 1.0 - term
        if abs(factor) < 1e-14:
            factor = complex(1e-14, 0.0)
        prod *= factor
    return 1.0 / prod


def evaluate_curve(
    cycles: list[CycleData],
    *,
    s_min: float,
    s_max: float,
    s_points: int,
    multiplier_floor: float,
    max_period: int | None,
    use_multiplier_abs: bool,
) -> tuple[np.ndarray, np.ndarray]:
    s_values = np.linspace(float(s_min), float(s_max), int(s_points), dtype=float)
    zeta = np.empty_like(s_values, dtype=np.complex128)
    for i, s in enumerate(s_values):
        zeta[i] = finite_product_zeta(
            complex(s, 0.0),
            cycles,
            multiplier_floor=multiplier_floor,
            max_period=max_period,
            use_multiplier_abs=use_multiplier_abs,
        )
    return s_values, zeta


def write_curve_csv(path: Path, s_values: np.ndarray, zeta: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["s", "zeta_real", "zeta_imag", "abs_zeta", "log10_abs_zeta"])
        for s, z in zip(s_values, zeta):
            abs_z = abs(z)
            w.writerow(
                [
                    f"{float(s):.16g}",
                    f"{z.real:.16g}",
                    f"{z.imag:.16g}",
                    f"{abs_z:.16g}",
                    f"{np.log10(max(abs_z, 1e-300)):.16g}",
                ]
            )


def plot_curve(path: Path, s_values: np.ndarray, zeta: np.ndarray, title: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(8.0, 4.8))
    plt.plot(s_values, np.abs(zeta), lw=1.8)
    plt.xlabel("s")
    plt.ylabel("|zeta(s)|")
    plt.title(title)
    plt.grid(alpha=0.25)
    plt.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Evaluate finite-product dynamical zeta approximation from cycle table."
    )
    parser.add_argument("--cycles", required=True, help="Path to cycle text table (from find_cycles.py).")
    parser.add_argument("--s-min", type=float, default=0.0)
    parser.add_argument("--s-max", type=float, default=3.0)
    parser.add_argument("--s-points", type=int, default=600)
    parser.add_argument("--multiplier-floor", type=float, default=1e-6)
    parser.add_argument("--max-period", type=int, default=0, help="0 means no cutoff.")
    parser.add_argument(
        "--use-multiplier-abs",
        action="store_true",
        help="Use |multiplier| in product terms instead of complex multiplier.",
    )
    parser.add_argument("--out-csv", default="outputs/csv/dynamical_zeta/zeta_curve.csv")
    parser.add_argument("--plot-out", default="outputs/figures/dynamical_zeta_abs_zeta.png")
    args = parser.parse_args()

    cycles_path = Path(args.cycles).resolve()
    cycles = read_cycle_table(cycles_path)
    if not cycles:
        raise RuntimeError(f"No cycles found in {cycles_path}")

    max_period = None if int(args.max_period) <= 0 else int(args.max_period)
    s_values, zeta = evaluate_curve(
        cycles,
        s_min=args.s_min,
        s_max=args.s_max,
        s_points=args.s_points,
        multiplier_floor=args.multiplier_floor,
        max_period=max_period,
        use_multiplier_abs=bool(args.use_multiplier_abs),
    )

    out_csv = Path(args.out_csv).resolve()
    plot_out = Path(args.plot_out).resolve()
    write_curve_csv(out_csv, s_values, zeta)
    plot_curve(plot_out, s_values, zeta, title=f"|zeta(s)| from cycles ({cycles_path.name})")

    print("Finite-product zeta")
    print(f"  cycles      : {cycles_path}")
    print(f"  cycle count : {len(cycles)}")
    print(f"  s range     : [{args.s_min}, {args.s_max}] with {args.s_points} points")
    print(f"  out csv     : {out_csv}")
    print(f"  plot        : {plot_out}")


if __name__ == "__main__":
    main()

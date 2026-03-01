from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# Allow running as: python scripts/dynamical_zeta/find_cycles.py
ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from collatzscape.maps import step_pair


@dataclass
class CycleRecord:
    period: int
    multiplier: complex
    location: complex
    hits: int = 1

    @property
    def multiplier_abs(self) -> float:
        try:
            val = float(abs(self.multiplier))
        except OverflowError:
            return float("inf")
        if not np.isfinite(val):
            return float("inf")
        return val

    @property
    def stability(self) -> str:
        m = self.multiplier_abs
        if not np.isfinite(m):
            return "unknown"
        if m < 1.0 - 1e-3:
            return "attracting"
        if m > 1.0 + 1e-3:
            return "repelling"
        return "neutral"


def make_grid(
    re_min: float,
    re_max: float,
    im_min: float,
    im_max: float,
    n_re: int,
    n_im: int,
) -> np.ndarray:
    x = np.linspace(re_min, re_max, int(n_re), dtype=float)
    y = np.linspace(im_min, im_max, int(n_im), dtype=float)
    gx, gy = np.meshgrid(x, y)
    return (gx + 1j * gy).ravel()


def _is_finite_complex(z: complex) -> bool:
    return np.isfinite(z.real) and np.isfinite(z.imag)


def _iterate_pair(
    z0: complex,
    *,
    steps: int,
    c: complex,
    alpha: float,
    gamma: float,
    direction: int,
    kick_mode: str,
    C_mode: str,
    escape_radius: float,
) -> tuple[list[complex], list[complex], bool]:
    z = complex(z0)
    w = complex(np.conjugate(z0))
    zs: list[complex] = []
    ws: list[complex] = []

    for _ in range(int(steps)):
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            z, w = step_pair(
                z,
                w,
                c=c,
                alpha=alpha,
                gamma=gamma,
                direction=direction,
                kick_mode=kick_mode,
                C_mode=C_mode,
            )

        if (not _is_finite_complex(z)) or (not _is_finite_complex(w)):
            return zs, ws, True
        if abs(z) > escape_radius or abs(w) > escape_radius:
            return zs, ws, True

        zs.append(z)
        ws.append(w)

    return zs, ws, False


def _detect_period(
    hist_z: list[complex],
    hist_w: list[complex],
    *,
    max_period: int,
    confirm_depth: int,
    close_tol: float,
) -> int | None:
    n = len(hist_z)
    if n < 3:
        return None

    for p in range(1, max_period + 1):
        if (p + confirm_depth) >= n:
            break
        ok = True
        for k in range(confirm_depth):
            if abs(hist_z[-1 - k] - hist_z[-1 - p - k]) > close_tol:
                ok = False
                break
            if abs(hist_w[-1 - k] - hist_w[-1 - p - k]) > close_tol:
                ok = False
                break
        if ok:
            return p
    return None


def _p_step_return_z(
    z0: complex,
    *,
    period: int,
    c: complex,
    alpha: float,
    gamma: float,
    direction: int,
    kick_mode: str,
    C_mode: str,
    escape_radius: float,
) -> complex:
    z = complex(z0)
    w = complex(np.conjugate(z0))
    for _ in range(int(period)):
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            z, w = step_pair(
                z,
                w,
                c=c,
                alpha=alpha,
                gamma=gamma,
                direction=direction,
                kick_mode=kick_mode,
                C_mode=C_mode,
            )
        if (not _is_finite_complex(z)) or abs(z) > escape_radius:
            return complex(np.nan, np.nan)
    return z


def estimate_cycle_multiplier(
    z_rep: complex,
    *,
    period: int,
    c: complex,
    alpha: float,
    gamma: float,
    direction: int,
    kick_mode: str,
    C_mode: str,
    escape_radius: float,
    eps: float = 1e-6,
) -> complex:
    """
    Multiplier estimate via dominant eigenvalue of finite-difference Jacobian
    of the p-step return map z -> F^p(z).
    """
    fr_p = _p_step_return_z(
        z_rep + eps,
        period=period,
        c=c,
        alpha=alpha,
        gamma=gamma,
        direction=direction,
        kick_mode=kick_mode,
        C_mode=C_mode,
        escape_radius=escape_radius,
    )
    fr_m = _p_step_return_z(
        z_rep - eps,
        period=period,
        c=c,
        alpha=alpha,
        gamma=gamma,
        direction=direction,
        kick_mode=kick_mode,
        C_mode=C_mode,
        escape_radius=escape_radius,
    )
    fi_p = _p_step_return_z(
        z_rep + 1j * eps,
        period=period,
        c=c,
        alpha=alpha,
        gamma=gamma,
        direction=direction,
        kick_mode=kick_mode,
        C_mode=C_mode,
        escape_radius=escape_radius,
    )
    fi_m = _p_step_return_z(
        z_rep - 1j * eps,
        period=period,
        c=c,
        alpha=alpha,
        gamma=gamma,
        direction=direction,
        kick_mode=kick_mode,
        C_mode=C_mode,
        escape_radius=escape_radius,
    )

    vals = [fr_p, fr_m, fi_p, fi_m]
    if any((not _is_finite_complex(v)) for v in vals):
        return complex(np.nan, np.nan)

    dzdx = (fr_p - fr_m) / (2.0 * eps)
    dzdy = (fi_p - fi_m) / (2.0 * eps)
    J = np.array([[dzdx.real, dzdy.real], [dzdx.imag, dzdy.imag]], dtype=float)
    if not np.all(np.isfinite(J)):
        return complex(np.nan, np.nan)

    eigvals = np.linalg.eigvals(J)
    idx = int(np.argmax(np.abs(eigvals)))
    lam = eigvals[idx]
    return complex(lam)


def discover_cycles(
    *,
    c: complex,
    alpha: float,
    gamma: float,
    direction: int,
    kick_mode: str,
    C_mode: str,
    re_min: float,
    re_max: float,
    im_min: float,
    im_max: float,
    n_re: int,
    n_im: int,
    steps: int,
    burn_in: int,
    max_period: int,
    confirm_depth: int,
    close_tol: float,
    dedup_tol: float,
    escape_radius: float,
    max_cycles: int = 10_000,
    verbose: bool = False,
) -> tuple[list[CycleRecord], dict[str, int]]:
    seeds = make_grid(re_min, re_max, im_min, im_max, n_re, n_im)
    cycles: list[CycleRecord] = []
    stats = {
        "seeds_total": int(seeds.size),
        "seeds_escaped": 0,
        "seeds_cycle_candidate": 0,
        "cycles_unique": 0,
    }

    max_hist = int(max_period + confirm_depth + 3)

    for idx, z0 in enumerate(seeds):
        z = complex(z0)
        w = complex(np.conjugate(z0))
        hist_z: list[complex] = []
        hist_w: list[complex] = []
        escaped = False
        detected_period: int | None = None

        for t in range(int(steps)):
            with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                z, w = step_pair(
                    z,
                    w,
                    c=c,
                    alpha=alpha,
                    gamma=gamma,
                    direction=direction,
                    kick_mode=kick_mode,
                    C_mode=C_mode,
                )

            if (not _is_finite_complex(z)) or (not _is_finite_complex(w)) or abs(z) > escape_radius or abs(w) > escape_radius:
                escaped = True
                break

            hist_z.append(z)
            hist_w.append(w)
            if len(hist_z) > max_hist:
                hist_z.pop(0)
                hist_w.pop(0)

            if t < burn_in:
                continue

            p = _detect_period(
                hist_z,
                hist_w,
                max_period=max_period,
                confirm_depth=confirm_depth,
                close_tol=close_tol,
            )
            if p is not None:
                detected_period = p
                break

        if escaped:
            stats["seeds_escaped"] += 1
            continue
        if detected_period is None:
            continue

        stats["seeds_cycle_candidate"] += 1
        p = int(detected_period)
        cycle_slice_z = np.array(hist_z[-p:], dtype=np.complex128)
        location = complex(np.mean(cycle_slice_z))
        z_rep = complex(cycle_slice_z[0])
        multiplier = estimate_cycle_multiplier(
            z_rep,
            period=p,
            c=c,
            alpha=alpha,
            gamma=gamma,
            direction=direction,
            kick_mode=kick_mode,
            C_mode=C_mode,
            escape_radius=escape_radius,
        )
        candidate = CycleRecord(period=p, multiplier=multiplier, location=location, hits=1)

        merged = False
        for rec in cycles:
            if rec.period != candidate.period:
                continue
            if abs(rec.location - candidate.location) <= dedup_tol:
                total_hits = rec.hits + 1
                rec.location = (rec.location * rec.hits + candidate.location) / total_hits
                if _is_finite_complex(rec.multiplier) and _is_finite_complex(candidate.multiplier):
                    rec.multiplier = (rec.multiplier * rec.hits + candidate.multiplier) / total_hits
                rec.hits = total_hits
                merged = True
                break
        if not merged:
            cycles.append(candidate)

        if verbose and (idx + 1) % 500 == 0:
            print(
                f"[find_cycles] processed {idx + 1}/{seeds.size}, "
                f"candidates={stats['seeds_cycle_candidate']}, unique={len(cycles)}"
            )
        if len(cycles) >= max_cycles:
            break

    cycles.sort(key=lambda r: (r.period, r.multiplier_abs, r.location.real, r.location.imag))
    stats["cycles_unique"] = len(cycles)
    return cycles, stats


def write_cycle_table(
    path: Path,
    cycles: list[CycleRecord],
    *,
    c: complex,
    alpha: float,
    gamma: float,
    direction: int,
    kick_mode: str,
    C_mode: str,
    stats: dict[str, int],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        f.write("# Dynamical cycle table for finite-product zeta approximation\n")
        f.write(f"# c={c}\n")
        f.write(f"# alpha={alpha}, gamma={gamma}, direction={direction}, kick_mode={kick_mode}, C_mode={C_mode}\n")
        f.write(
            f"# seeds_total={stats.get('seeds_total', 0)}, seeds_escaped={stats.get('seeds_escaped', 0)}, "
            f"cycle_candidates={stats.get('seeds_cycle_candidate', 0)}, cycles_unique={stats.get('cycles_unique', 0)}\n"
        )

        writer = csv.writer(f, delimiter="\t")
        writer.writerow(
            [
                "id",
                "period",
                "multiplier_real",
                "multiplier_imag",
                "multiplier_abs",
                "location_real",
                "location_imag",
                "hits",
                "stability",
            ]
        )
        for i, rec in enumerate(cycles):
            writer.writerow(
                [
                    i,
                    rec.period,
                    f"{rec.multiplier.real:.16g}",
                    f"{rec.multiplier.imag:.16g}",
                    f"{rec.multiplier_abs:.16g}",
                    f"{rec.location.real:.16g}",
                    f"{rec.location.imag:.16g}",
                    rec.hits,
                    rec.stability,
                ]
            )


def parse_complex(raw: str) -> complex:
    return complex(raw.strip().replace("i", "j"))


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate grid seeds, iterate the Collatz-tetration pair map, detect close returns, "
            "and write cycle data (period, multiplier, location) to text."
        )
    )
    parser.add_argument("--c", default="1.93+0j")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--gamma", type=float, default=1.0)
    parser.add_argument("--direction", type=int, choices=[-1, 1], default=1)
    parser.add_argument("--kick-mode", choices=["tanh", "clip"], default="tanh")
    parser.add_argument("--C-mode", choices=["cosine", "linear"], default="cosine")

    parser.add_argument("--re-min", type=float, default=-2.0)
    parser.add_argument("--re-max", type=float, default=2.0)
    parser.add_argument("--im-min", type=float, default=-2.0)
    parser.add_argument("--im-max", type=float, default=2.0)
    parser.add_argument("--n-re", type=int, default=110)
    parser.add_argument("--n-im", type=int, default=110)

    parser.add_argument("--steps", type=int, default=220)
    parser.add_argument("--burn-in", type=int, default=80)
    parser.add_argument("--max-period", type=int, default=24)
    parser.add_argument("--confirm-depth", type=int, default=4)
    parser.add_argument("--close-tol", type=float, default=1e-5)
    parser.add_argument("--dedup-tol", type=float, default=3e-4)
    parser.add_argument("--escape-radius", type=float, default=1e6)
    parser.add_argument("--max-cycles", type=int, default=2000)

    parser.add_argument("--out", default="outputs/csv/dynamical_zeta/cycles.txt")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    c = parse_complex(args.c)
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

    out = Path(args.out).resolve()
    write_cycle_table(
        out,
        cycles,
        c=c,
        alpha=args.alpha,
        gamma=args.gamma,
        direction=args.direction,
        kick_mode=args.kick_mode,
        C_mode=args.C_mode,
        stats=stats,
    )

    print("Dynamical cycle discovery")
    print(f"  c                 : {c}")
    print(f"  seeds             : {stats['seeds_total']}")
    print(f"  escaped seeds     : {stats['seeds_escaped']}")
    print(f"  cycle candidates  : {stats['seeds_cycle_candidate']}")
    print(f"  unique cycles     : {stats['cycles_unique']}")
    print(f"  output            : {out}")


if __name__ == "__main__":
    main()

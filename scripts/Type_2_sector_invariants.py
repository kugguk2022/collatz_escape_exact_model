from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

# Allow running as: python scripts/Type_2_sector_invariants.py
ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from collatzscape.maps import step_pair


TOWER_A_MIN = float(np.exp(-np.e)) + 1e-6
TOWER_A_MAX = float(np.exp(1.0 / np.e)) - 1e-6
ABS_FLOOR = 1e-15


@dataclass
class StepRecord:
    t: int
    a_t: float
    x_tower: float
    kappa_t: float
    gain_t: float
    theta_t: float
    gamma_t: float
    gamma_next: float
    P_t: int
    M_t: int
    direction_t: int
    z_abs: float
    z_real: float
    z_imag: float
    escaped: int


@dataclass
class RunSummary:
    steps_requested: int
    steps_executed: int
    escaped: bool
    escape_step: int
    max_abs_z: float
    final_abs_z: float
    mean_abs_z: float
    lyapunov_proxy: float
    mean_gamma: float
    final_gamma: float
    mean_kappa: float
    mean_x_tower: float
    final_a: float
    final_P: int
    final_M: int


def clamp(value: float, lo: float, hi: float) -> float:
    return float(min(max(value, lo), hi))


def safe_abs(z: complex) -> float:
    if (not np.isfinite(z.real)) or (not np.isfinite(z.imag)):
        return math.inf
    return float(abs(z))


def sieve_primes_upto(n: int) -> list[int]:
    n = int(n)
    if n < 2:
        return []
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    limit = int(np.sqrt(n))
    for i in range(2, limit + 1):
        if s[i]:
            s[i * i : n + 1 : i] = False
    return np.flatnonzero(s).astype(int).tolist()


def bernoulli_zp(p: int, beta: float) -> float:
    """
    Bernoulli-style Z(p) map with deterministic convergence to 2:
        z_p = 2 - 1/(p+1)^beta
    """
    return float(2.0 - 1.0 / ((p + 1.0) ** beta))


def solve_power_tower_fixed_point(
    a_t: float,
    *,
    x0: float | None = None,
    max_iter: int = 120,
    tol: float = 1e-12,
    damping: float = 0.65,
) -> tuple[float, bool, int]:
    """
    Solve x = a^x by damped fixed-point iteration + Newton polish.
    Stable in the tower interval a in [exp(-e), exp(1/e)].
    """
    a = float(a_t)
    if (a <= 0.0) or (not np.isfinite(a)):
        return float("nan"), False, 0

    x = float(x0 if x0 is not None else max(1.0, a))
    converged = False

    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        for i in range(max_iter):
            x_raw = float(a**x)
            if not np.isfinite(x_raw):
                return float("nan"), False, i
            x_next = (1.0 - damping) * x + damping * x_raw
            if not np.isfinite(x_next):
                return float("nan"), False, i
            if abs(x_next - x) <= tol * (1.0 + abs(x_next)):
                x = x_next
                converged = True
                break
            x = x_next

        # One short Newton polish pass helps numerical stability.
        for _ in range(10):
            ax = float(a**x)
            f = x - ax
            fp = 1.0 - math.log(a) * ax
            if abs(fp) < 1e-14:
                break
            x_new = x - f / fp
            if not np.isfinite(x_new):
                break
            if abs(x_new - x) <= tol * (1.0 + abs(x_new)):
                x = x_new
                converged = True
                break
            x = x_new

    return float(x), converged, i + 1


def prime_sector_kappa(
    x_tower: float,
    p_cap: int,
    *,
    beta_zp: float,
    regularization_scale: float,
) -> tuple[float, int]:
    """
    kappa = Reg( sum_{p<=P} w(p;x) * kernel(p;x) )

    Uses:
      w(p;x)      = z_p * p^{-x}
      z_p         = 2 - 1/(p+1)^beta_zp
      kernel(p;x) = -log(1 - p^{-x})    (Euler-factor style)
      Reg(y)      = tanh(regularization_scale * y)
    """
    primes = sieve_primes_upto(int(p_cap))
    if not primes:
        return 0.0, 0

    x = max(1.0000001, float(x_tower))
    total = 0.0

    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        for p in primes:
            px = float(p ** (-x))
            zp = bernoulli_zp(p, beta=beta_zp)
            weight = zp * px
            kernel = float(-math.log(max(1.0 - px, ABS_FLOOR)))
            total += weight * kernel

    kappa = float(np.tanh(regularization_scale * total))
    return kappa, len(primes)


def map_gain_from_kappa(kappa: float, gain0: float, gain_kappa: float, *, a_min: float, a_max: float) -> float:
    return clamp(gain0 + gain_kappa * math.tanh(kappa), a_min, a_max)


def map_theta_from_gamma(gamma_t: float, theta0: float, theta_gamma: float) -> float:
    return float(theta0 + theta_gamma * math.atan(gamma_t))


def catalan_numbers(n: int) -> np.ndarray:
    """
    C_0=1, C_{n+1} = C_n * 2(2n+1)/(n+2)
    """
    n = max(1, int(n))
    out = np.empty(n, dtype=float)
    out[0] = 1.0
    for k in range(1, n):
        out[k] = out[k - 1] * (2.0 * (2.0 * (k - 1) + 1.0)) / (k + 1.0)
    return out


def build_signal_from_z(z: complex, m_t: int) -> np.ndarray:
    """
    f(theta) = sign(Im(z * exp(i theta))) sampled over theta-grid.
    """
    m = max(8, int(m_t))
    theta = np.linspace(0.0, 2.0 * np.pi, m, endpoint=False)
    vals = np.imag(z * np.exp(1j * theta))
    f = np.sign(vals)
    f[f == 0] = 1.0
    return f.astype(float)


def odd_harmonic_energy(signal: np.ndarray, *, beta_catalan: float) -> float:
    """
    Odd-sector energy using Catalan-weighted odd Fourier modes.

    beta_catalan is constrained to >= 2 by caller.
    """
    m = int(signal.size)
    coeff = np.fft.rfft(signal) / float(m)
    odd_idx = np.arange(1, coeff.size, 2)
    if odd_idx.size == 0:
        return 0.0

    cats = catalan_numbers(odd_idx.size)
    weights = 1.0 / np.power(cats + 1.0, beta_catalan / 2.0)
    energy = float(np.sum(weights * (np.abs(coeff[odd_idx]) ** 2)))
    return max(0.0, energy)


def collatz_tetration_step(
    z_t: complex,
    w_t: complex,
    *,
    gain_t: float,
    theta_t: float,
    gamma_t: float,
    alpha: float,
    gamma_scale: float,
    gamma_clip: float,
    kick_mode: str,
    C_mode: str,
) -> tuple[complex, complex, int, float]:
    """
    Adapter from (gain,theta,gamma) to repo-native step_pair API.
    """
    c_t = complex(gain_t * np.exp(1j * theta_t))
    direction_t = +1 if theta_t >= 0.0 else -1
    gamma_model = clamp(gamma_scale * abs(gamma_t), 0.0, gamma_clip)

    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        z_next, w_next = step_pair(
            z_t,
            w_t,
            c=c_t,
            alpha=alpha,
            gamma=gamma_model,
            direction=direction_t,
            kick_mode=kick_mode,
            C_mode=C_mode,
        )
    return z_next, w_next, direction_t, gamma_model


def run_type2_sector_invariants(
    *,
    steps: int,
    z0: complex,
    a0: float,
    P0: int,
    M0: int,
    gain0: float,
    theta0: float,
    gamma0: float,
    A: float,
    B: float,
    gain_kappa: float,
    theta_gamma: float,
    alpha: float,
    beta_zp: float,
    beta_catalan: float,
    regularization_scale: float,
    gamma_scale: float,
    gamma_clip: float,
    escape_radius: float,
    kick_mode: str,
    C_mode: str,
    tower_max_iter: int,
    tower_tol: float,
    tower_damping: float,
    a_min: float,
    a_max: float,
) -> tuple[list[StepRecord], RunSummary]:
    if beta_catalan < 2.0:
        raise ValueError("beta_catalan must be >= 2.0")

    z_t = complex(z0)
    w_t = complex(np.conjugate(z0))
    a_t = clamp(float(a0), a_min, a_max)
    P_t = max(2, int(P0))
    M_t = max(8, int(M0))
    gamma_t = float(gamma0)

    records: list[StepRecord] = []
    growth_logs: list[float] = []
    kappa_vals: list[float] = []
    x_vals: list[float] = []
    gamma_vals: list[float] = []

    escaped = False
    escape_step = -1
    max_abs_z = safe_abs(z_t)

    for t in range(int(steps)):
        x_tower, _ok, _iters = solve_power_tower_fixed_point(
            a_t,
            x0=max(1.0, a_t),
            max_iter=tower_max_iter,
            tol=tower_tol,
            damping=tower_damping,
        )
        if not np.isfinite(x_tower):
            x_tower = max(1.0, a_t)

        kappa_t, _n_primes = prime_sector_kappa(
            x_tower,
            P_t,
            beta_zp=beta_zp,
            regularization_scale=regularization_scale,
        )

        gain_t = map_gain_from_kappa(kappa_t, gain0, gain_kappa, a_min=a_min, a_max=a_max)
        theta_t = map_theta_from_gamma(gamma_t, theta0, theta_gamma)

        z_prev_abs = safe_abs(z_t)
        z_t, w_t, direction_t, _gamma_model = collatz_tetration_step(
            z_t,
            w_t,
            gain_t=gain_t,
            theta_t=theta_t,
            gamma_t=gamma_t,
            alpha=alpha,
            gamma_scale=gamma_scale,
            gamma_clip=gamma_clip,
            kick_mode=kick_mode,
            C_mode=C_mode,
        )
        z_abs = safe_abs(z_t)
        max_abs_z = max(max_abs_z, z_abs)

        signal_t = build_signal_from_z(z_t, M_t)
        gamma_next = odd_harmonic_energy(signal_t, beta_catalan=beta_catalan)

        records.append(
            StepRecord(
                t=t,
                a_t=float(a_t),
                x_tower=float(x_tower),
                kappa_t=float(kappa_t),
                gain_t=float(gain_t),
                theta_t=float(theta_t),
                gamma_t=float(gamma_t),
                gamma_next=float(gamma_next),
                P_t=int(P_t),
                M_t=int(M_t),
                direction_t=int(direction_t),
                z_abs=float(z_abs),
                z_real=float(np.real(z_t)) if np.isfinite(z_t.real) else float("nan"),
                z_imag=float(np.imag(z_t)) if np.isfinite(z_t.imag) else float("nan"),
                escaped=int((not np.isfinite(z_abs)) or (z_abs > escape_radius)),
            )
        )

        kappa_vals.append(float(kappa_t))
        x_vals.append(float(x_tower))
        gamma_vals.append(float(gamma_next))

        if np.isfinite(z_prev_abs) and np.isfinite(z_abs):
            growth_logs.append(float(math.log((z_abs + ABS_FLOOR) / (z_prev_abs + ABS_FLOOR))))

        if (not np.isfinite(z_abs)) or (z_abs > escape_radius):
            escaped = True
            escape_step = t
            break

        # Renormalize truncations
        P_t = max(2, int(P0 + math.floor(A * x_tower)))
        M_t = max(8, int(M0 + math.floor(B * x_tower)))

        # Update latent states for next loop
        gamma_t = float(gamma_next)
        a_t = clamp(gain_t, a_min, a_max)

    steps_executed = len(records)
    mean_abs_z = float(np.mean([r.z_abs for r in records])) if records else float("nan")
    lyapunov_proxy = float(np.mean(growth_logs)) if growth_logs else float("nan")
    mean_gamma = float(np.mean(gamma_vals)) if gamma_vals else float("nan")
    mean_kappa = float(np.mean(kappa_vals)) if kappa_vals else float("nan")
    mean_x_tower = float(np.mean(x_vals)) if x_vals else float("nan")

    summary = RunSummary(
        steps_requested=int(steps),
        steps_executed=steps_executed,
        escaped=escaped,
        escape_step=int(escape_step),
        max_abs_z=float(max_abs_z),
        final_abs_z=float(records[-1].z_abs if records else safe_abs(z_t)),
        mean_abs_z=mean_abs_z,
        lyapunov_proxy=lyapunov_proxy,
        mean_gamma=mean_gamma,
        final_gamma=float(records[-1].gamma_next if records else gamma_t),
        mean_kappa=mean_kappa,
        mean_x_tower=mean_x_tower,
        final_a=float(a_t),
        final_P=int(P_t),
        final_M=int(M_t),
    )
    return records, summary


def write_csv(path: Path, records: list[StepRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(records[0]).keys()) if records else [])
        if records:
            w.writeheader()
            for rec in records:
                w.writerow(asdict(rec))


def write_summary_json(path: Path, summary: RunSummary) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(asdict(summary), f, indent=2)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Type-2 sector invariants: tower fixed point + prime sector + "
            "Collatz tetration + Catalan odd-sector energy."
        )
    )
    parser.add_argument("--steps", type=int, default=80)

    parser.add_argument("--z0-real", type=float, default=0.7)
    parser.add_argument("--z0-imag", type=float, default=0.2)
    parser.add_argument("--a0", type=float, default=float(np.sqrt(2.0)))
    parser.add_argument("--P0", type=int, default=64)
    parser.add_argument("--M0", type=int, default=64)

    parser.add_argument("--gain0", type=float, default=float(np.sqrt(2.0)))
    parser.add_argument("--theta0", type=float, default=0.0)
    parser.add_argument("--gamma0", type=float, default=0.0)

    parser.add_argument("--A", type=float, default=24.0)
    parser.add_argument("--B", type=float, default=24.0)
    parser.add_argument("--gain-kappa", type=float, default=0.08)
    parser.add_argument("--theta-gamma", type=float, default=0.9)
    parser.add_argument("--alpha", type=float, default=0.02)
    parser.add_argument("--gamma-scale", type=float, default=1.0)
    parser.add_argument("--gamma-clip", type=float, default=2.0)

    parser.add_argument("--beta-zp", type=float, default=2.0, help="Bernoulli Z(p) exponent.")
    parser.add_argument(
        "--beta-catalan",
        type=float,
        default=2.0,
        help="Catalan odd-sector exponent (must be >= 2).",
    )
    parser.add_argument("--regularization-scale", type=float, default=3.0)

    parser.add_argument("--kick-mode", choices=["tanh", "clip"], default="tanh")
    parser.add_argument("--C-mode", choices=["cosine", "linear"], default="cosine")
    parser.add_argument("--escape-radius", type=float, default=1e6)

    parser.add_argument("--tower-max-iter", type=int, default=120)
    parser.add_argument("--tower-tol", type=float, default=1e-12)
    parser.add_argument("--tower-damping", type=float, default=0.65)
    parser.add_argument("--a-min", type=float, default=TOWER_A_MIN)
    parser.add_argument("--a-max", type=float, default=TOWER_A_MAX)

    parser.add_argument("--out-csv", default="outputs/csv/type_2_sector_invariants.csv")
    parser.add_argument("--summary-json", default="outputs/csv/type_2_sector_invariants_summary.json")
    args = parser.parse_args()

    records, summary = run_type2_sector_invariants(
        steps=args.steps,
        z0=complex(args.z0_real, args.z0_imag),
        a0=args.a0,
        P0=args.P0,
        M0=args.M0,
        gain0=args.gain0,
        theta0=args.theta0,
        gamma0=args.gamma0,
        A=args.A,
        B=args.B,
        gain_kappa=args.gain_kappa,
        theta_gamma=args.theta_gamma,
        alpha=args.alpha,
        beta_zp=args.beta_zp,
        beta_catalan=args.beta_catalan,
        regularization_scale=args.regularization_scale,
        gamma_scale=args.gamma_scale,
        gamma_clip=args.gamma_clip,
        escape_radius=args.escape_radius,
        kick_mode=args.kick_mode,
        C_mode=args.C_mode,
        tower_max_iter=args.tower_max_iter,
        tower_tol=args.tower_tol,
        tower_damping=args.tower_damping,
        a_min=args.a_min,
        a_max=args.a_max,
    )

    out_csv = Path(args.out_csv).resolve()
    out_summary = Path(args.summary_json).resolve()
    if records:
        write_csv(out_csv, records)
    write_summary_json(out_summary, summary)

    print("Type-2 sector invariants")
    print(f"  steps_executed : {summary.steps_executed}/{summary.steps_requested}")
    print(f"  escaped        : {summary.escaped} at step {summary.escape_step}")
    print(f"  max_abs_z      : {summary.max_abs_z:.6g}")
    print(f"  final_abs_z    : {summary.final_abs_z:.6g}")
    print(f"  lyapunov_proxy : {summary.lyapunov_proxy:.6g}")
    print(f"  final_gamma    : {summary.final_gamma:.6g}")
    print(f"  mean_x_tower   : {summary.mean_x_tower:.6g}")
    print(f"  final_P, final_M: {summary.final_P}, {summary.final_M}")
    print(f"  csv            : {out_csv}")
    print(f"  summary_json   : {out_summary}")


if __name__ == "__main__":
    main()

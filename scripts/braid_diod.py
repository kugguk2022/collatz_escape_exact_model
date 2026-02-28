from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# Allow running as: python scripts/braid_diod.py
ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from collatzscape.maps import C_cosine, C_linear_near2, c_eff, f_clip, f_tanh, step_pair


LOG_FLOOR = 1e-15


def _safe_abs(z: complex) -> float:
    if (not np.isfinite(z.real)) or (not np.isfinite(z.imag)):
        return math.inf
    return float(abs(z))


def _gate_signal(z: complex, *, kick_mode: str, C_mode: str) -> float:
    c_fun = C_cosine if C_mode == "cosine" else C_linear_near2
    x = float(np.real(c_fun(z)))
    if kick_mode == "tanh":
        return f_tanh(x)
    if kick_mode == "clip":
        return f_clip(x)
    raise ValueError(f"Unknown kick mode: {kick_mode}")


@dataclass
class DiodeObservation:
    delta_f: float
    log_ratio: float
    abs_c_eff_plus: float
    abs_c_eff_minus: float
    valid: bool


def observe_diode(
    z: complex,
    w: complex,
    *,
    c: complex,
    alpha: float,
    gamma_true: float,
    kick_mode: str,
    C_mode: str,
) -> DiodeObservation:
    """
    Exact diode observable from current (z,w):

        y = log|c_eff(dir=-1)| - log|c_eff(dir=+1)| = gamma * delta_f

    where delta_f = f(Re(C(z))) - f(Re(C(w))).
    """
    sig_z = _gate_signal(z, kick_mode=kick_mode, C_mode=C_mode)
    sig_w = _gate_signal(w, kick_mode=kick_mode, C_mode=C_mode)
    delta_f = float(sig_z - sig_w)

    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        ce_plus = c_eff(
            z,
            w,
            c=c,
            alpha=alpha,
            gamma=gamma_true,
            direction=+1,
            kick_mode=kick_mode,
            C_mode=C_mode,
        )
        ce_minus = c_eff(
            z,
            w,
            c=c,
            alpha=alpha,
            gamma=gamma_true,
            direction=-1,
            kick_mode=kick_mode,
            C_mode=C_mode,
        )

    abs_plus = _safe_abs(ce_plus)
    abs_minus = _safe_abs(ce_minus)

    if np.isfinite(abs_plus) and np.isfinite(abs_minus):
        log_ratio = float(math.log(max(abs_minus, LOG_FLOOR)) - math.log(max(abs_plus, LOG_FLOOR)))
        valid = np.isfinite(log_ratio) and np.isfinite(delta_f)
    else:
        log_ratio = float("nan")
        valid = False

    return DiodeObservation(
        delta_f=delta_f,
        log_ratio=log_ratio,
        abs_c_eff_plus=abs_plus,
        abs_c_eff_minus=abs_minus,
        valid=valid,
    )


@dataclass
class OnlineGammaEstimator:
    reg_eps: float = 1e-12
    s_xx: float = 0.0
    s_xy: float = 0.0
    gamma_hat: float = 0.0
    samples_used: int = 0

    def update(self, x: float, y: float) -> float:
        if (not np.isfinite(x)) or (not np.isfinite(y)):
            return self.gamma_hat
        self.s_xx += float(x * x)
        self.s_xy += float(x * y)
        if self.s_xx > self.reg_eps:
            self.gamma_hat = float(self.s_xy / self.s_xx)
            self.samples_used += 1
        return self.gamma_hat


@dataclass
class RunSummary:
    mode: str
    escaped: bool
    escape_step: int
    steps_executed: int
    max_abs_z: float
    final_abs_z: float
    gamma_hat: float
    gamma_abs_error: float
    controller_plus_fraction: float
    controller_switches: int
    estimator_samples: int


def choose_direction_lookahead(
    z: complex,
    w: complex,
    *,
    c: complex,
    alpha: float,
    gamma_for_prediction: float,
    kick_mode: str,
    C_mode: str,
    fallback: int,
) -> int:
    """
    One-step stabilizing controller:
    evaluate both directional updates and pick the branch with smaller |z_next|.
    """
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        zp, _ = step_pair(
            z,
            w,
            c=c,
            alpha=alpha,
            gamma=gamma_for_prediction,
            direction=+1,
            kick_mode=kick_mode,
            C_mode=C_mode,
        )
        zm, _ = step_pair(
            z,
            w,
            c=c,
            alpha=alpha,
            gamma=gamma_for_prediction,
            direction=-1,
            kick_mode=kick_mode,
            C_mode=C_mode,
        )

    ap = _safe_abs(zp)
    am = _safe_abs(zm)
    if (not np.isfinite(ap)) and (not np.isfinite(am)):
        return int(fallback)
    return +1 if ap <= am else -1


def run_episode(
    *,
    mode: str,
    z0: complex,
    c: complex,
    alpha: float,
    gamma_true: float,
    steps: int,
    escape_radius: float,
    open_loop_direction: int,
    kick_mode: str,
    C_mode: str,
    gamma_pred_floor: float,
    trace_rows: list[dict[str, float | int | str]] | None = None,
) -> RunSummary:
    estimator = OnlineGammaEstimator()

    z = complex(z0)
    w = complex(np.conjugate(z0))
    max_abs_z = _safe_abs(z)
    escaped = False
    escape_step = -1

    last_direction = int(open_loop_direction)
    plus_count = 0
    n_switches = 0

    for t in range(steps):
        obs = observe_diode(
            z,
            w,
            c=c,
            alpha=alpha,
            gamma_true=gamma_true,
            kick_mode=kick_mode,
            C_mode=C_mode,
        )

        if obs.valid:
            gamma_hat = estimator.update(obs.delta_f, obs.log_ratio)
        else:
            gamma_hat = estimator.gamma_hat

        if mode == "closed_loop":
            gamma_pred = max(abs(gamma_hat), gamma_pred_floor)
            direction = choose_direction_lookahead(
                z,
                w,
                c=c,
                alpha=alpha,
                gamma_for_prediction=gamma_pred,
                kick_mode=kick_mode,
                C_mode=C_mode,
                fallback=last_direction,
            )
        else:
            direction = int(open_loop_direction)

        if direction == +1:
            plus_count += 1
        if t > 0 and direction != last_direction:
            n_switches += 1
        last_direction = direction

        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            z, w = step_pair(
                z,
                w,
                c=c,
                alpha=alpha,
                gamma=gamma_true,
                direction=direction,
                kick_mode=kick_mode,
                C_mode=C_mode,
            )

        abs_z = _safe_abs(z)
        max_abs_z = max(max_abs_z, abs_z)

        if trace_rows is not None:
            trace_rows.append(
                {
                    "mode": mode,
                    "step": t,
                    "direction": direction,
                    "abs_z": abs_z,
                    "abs_c_eff_plus": obs.abs_c_eff_plus,
                    "abs_c_eff_minus": obs.abs_c_eff_minus,
                    "delta_f": obs.delta_f,
                    "log_ratio": obs.log_ratio,
                    "gamma_hat": gamma_hat,
                }
            )

        if (not np.isfinite(abs_z)) or (abs_z > escape_radius):
            escaped = True
            escape_step = t
            break

    steps_executed = (escape_step + 1) if escaped else steps
    plus_fraction = (plus_count / max(1, steps_executed))
    gamma_hat_final = estimator.gamma_hat

    return RunSummary(
        mode=mode,
        escaped=escaped,
        escape_step=escape_step,
        steps_executed=steps_executed,
        max_abs_z=max_abs_z,
        final_abs_z=_safe_abs(z),
        gamma_hat=gamma_hat_final,
        gamma_abs_error=abs(gamma_hat_final - gamma_true),
        controller_plus_fraction=plus_fraction,
        controller_switches=n_switches,
        estimator_samples=estimator.samples_used,
    )


def write_trace_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "mode",
                "step",
                "direction",
                "abs_z",
                "abs_c_eff_plus",
                "abs_c_eff_minus",
                "delta_f",
                "log_ratio",
                "gamma_hat",
            ],
        )
        w.writeheader()
        w.writerows(rows)


def print_summary(title: str, summary: RunSummary) -> None:
    print(f"\n{title}")
    print(f"  mode                 : {summary.mode}")
    print(f"  escaped              : {summary.escaped}")
    print(f"  escape_step          : {summary.escape_step}")
    print(f"  steps_executed       : {summary.steps_executed}")
    print(f"  max_abs_z            : {summary.max_abs_z:.6g}")
    print(f"  final_abs_z          : {summary.final_abs_z:.6g}")
    print(f"  gamma_hat            : {summary.gamma_hat:.6g}")
    print(f"  gamma_abs_error      : {summary.gamma_abs_error:.6g}")
    print(f"  controller_plus_frac : {summary.controller_plus_fraction:.6g}")
    print(f"  controller_switches  : {summary.controller_switches}")
    print(f"  estimator_samples    : {summary.estimator_samples}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Braid-Diode: estimate gamma from exact forward/reverse asymmetry and "
            "apply a stabilizing direction controller."
        )
    )
    parser.add_argument("--c-real", type=float, default=1.5)
    parser.add_argument("--c-imag", type=float, default=0.5)
    parser.add_argument("--z0-real", type=float, default=0.7)
    parser.add_argument("--z0-imag", type=float, default=0.2)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--gamma", type=float, default=0.3, help="True gamma used by the simulated plant.")
    parser.add_argument("--steps", type=int, default=80)
    parser.add_argument("--escape-radius", type=float, default=1e6)
    parser.add_argument("--open-loop-direction", type=int, choices=[-1, 1], default=-1)
    parser.add_argument("--kick-mode", choices=["tanh", "clip"], default="tanh")
    parser.add_argument("--C-mode", choices=["cosine", "linear"], default="cosine")
    parser.add_argument(
        "--gamma-pred-floor",
        type=float,
        default=1e-6,
        help="Minimum |gamma| used for one-step look-ahead when gamma_hat is near zero.",
    )
    parser.add_argument(
        "--trace-csv",
        default="",
        help="Optional path for per-step diagnostics CSV (both open_loop and closed_loop rows).",
    )
    args = parser.parse_args()

    c = complex(args.c_real, args.c_imag)
    z0 = complex(args.z0_real, args.z0_imag)

    trace_rows: list[dict[str, float | int | str]] | None = [] if args.trace_csv else None

    open_loop = run_episode(
        mode="open_loop",
        z0=z0,
        c=c,
        alpha=args.alpha,
        gamma_true=args.gamma,
        steps=args.steps,
        escape_radius=args.escape_radius,
        open_loop_direction=args.open_loop_direction,
        kick_mode=args.kick_mode,
        C_mode=args.C_mode,
        gamma_pred_floor=args.gamma_pred_floor,
        trace_rows=trace_rows,
    )
    closed_loop = run_episode(
        mode="closed_loop",
        z0=z0,
        c=c,
        alpha=args.alpha,
        gamma_true=args.gamma,
        steps=args.steps,
        escape_radius=args.escape_radius,
        open_loop_direction=args.open_loop_direction,
        kick_mode=args.kick_mode,
        C_mode=args.C_mode,
        gamma_pred_floor=args.gamma_pred_floor,
        trace_rows=trace_rows,
    )

    print("Braid-Diode")
    print(f"  c                    : {c}")
    print(f"  z0                   : {z0}")
    print(f"  alpha                : {args.alpha}")
    print(f"  gamma_true           : {args.gamma}")
    print_summary("Open-loop baseline", open_loop)
    print_summary("Closed-loop (stabilizing)", closed_loop)

    if args.trace_csv:
        trace_path = Path(args.trace_csv).resolve()
        assert trace_rows is not None
        write_trace_csv(trace_path, trace_rows)
        print(f"\ntrace_csv: {trace_path}")


if __name__ == "__main__":
    main()

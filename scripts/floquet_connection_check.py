from __future__ import annotations

import argparse
import importlib.util
import math
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from collatzscape.maps import C_cosine, c_eff, f_tanh, step_pair


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module from {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def run_stroboscopic_and_gain_loss_check(
    *,
    c: complex,
    alpha: float,
    gamma: float,
    direction: int,
    z: complex,
    w: complex,
    kick_mode: str,
    C_mode: str,
) -> dict[str, object]:
    z1, w1 = step_pair(
        z,
        w,
        c=c,
        alpha=alpha,
        gamma=gamma,
        direction=direction,
        kick_mode=kick_mode,
        C_mode=C_mode,
    )

    sig_z = f_tanh(float(np.real(C_cosine(z)))) if kick_mode == "tanh" else float(np.clip(float(np.real(C_cosine(z))), -1.0, 1.0))
    sig_w = f_tanh(float(np.real(C_cosine(w)))) if kick_mode == "tanh" else float(np.clip(float(np.real(C_cosine(w))), -1.0, 1.0))
    delta_f = float(sig_z - sig_w)

    ce_plus = c_eff(
        z,
        w,
        c=c,
        alpha=alpha,
        gamma=gamma,
        direction=+1,
        kick_mode=kick_mode,
        C_mode=C_mode,
    )
    ce_minus = c_eff(
        z,
        w,
        c=c,
        alpha=alpha,
        gamma=gamma,
        direction=-1,
        kick_mode=kick_mode,
        C_mode=C_mode,
    )
    log_ratio = float(math.log(abs(ce_minus)) - math.log(abs(ce_plus)))
    identity_rhs = float(gamma * delta_f)

    return {
        "z0": z,
        "w0": w,
        "z1": z1,
        "w1": w1,
        "delta_f": delta_f,
        "log_ratio": log_ratio,
        "gamma_delta_f": identity_rhs,
        "identity_abs_error": float(abs(log_ratio - identity_rhs)),
    }


def run_cycle_multiplier_check(
    *,
    c: complex,
    alpha: float,
    gamma: float,
    direction: int,
    kick_mode: str,
    C_mode: str,
    out_path: Path,
) -> dict[str, object]:
    find_cycles = _load_module(ROOT / "scripts" / "dynamical_zeta" / "find_cycles.py", "find_cycles_mod")
    cycles, stats = find_cycles.discover_cycles(
        c=c,
        alpha=alpha,
        gamma=gamma,
        direction=direction,
        kick_mode=kick_mode,
        C_mode=C_mode,
        re_min=-2.0,
        re_max=2.0,
        im_min=-2.0,
        im_max=2.0,
        n_re=28,
        n_im=28,
        steps=120,
        burn_in=24,
        max_period=14,
        confirm_depth=3,
        close_tol=1e-5,
        dedup_tol=3e-4,
        escape_radius=1e6,
        max_cycles=6,
        verbose=False,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    find_cycles.write_cycle_table(
        out_path,
        cycles,
        c=c,
        alpha=alpha,
        gamma=gamma,
        direction=direction,
        kick_mode=kick_mode,
        C_mode=C_mode,
        stats=stats,
    )
    first = cycles[0] if cycles else None
    return {
        "stats": stats,
        "cycles_unique": int(len(cycles)),
        "first_period": int(first.period) if first else None,
        "first_multiplier_abs": float(first.multiplier_abs) if first else None,
        "first_stability": str(first.stability) if first else None,
        "out_file": str(out_path),
    }


def run_quasiparticle_mode_check(*, alpha: float) -> dict[str, float]:
    mod = _load_module(ROOT / "scripts" / "Collatz_tetration_quasiparticles.py", "quasi_mod")
    c_values = [1.8725 + 0j, 2.2 + 0j, 1.87 + 0.35j, -1.2 + 0.8j]
    out: dict[str, float] = {}
    for c in c_values:
        s = mod.survival_probability(
            c,
            z0=c,
            alpha=alpha,
            n_steps=120,
            radius=80.0,
            trials=24,
            dz=0.03,
            da=0.01,
            seed=7,
            phase_mode="tanh",
            tanh_scale=1.0,
        )
        out[str(c)] = float(s)
    return out


def build_markdown_report(
    stroboscopic: dict[str, object],
    cycle: dict[str, object],
    modes: dict[str, float],
) -> str:
    rows = [
        (
            "Stroboscopic update over one drive period",
            "one-step map `(z_n,w_n) -> (z_{n+1},w_{n+1})`",
            "`src/collatzscape/maps.py` (`step_pair`)",
            f"`z0={stroboscopic['z0']}, w0={stroboscopic['w0']} -> z1={stroboscopic['z1']}, w1={stroboscopic['w1']}`",
        ),
        (
            "Period-`p` orbit stability via multipliers",
            "cycle multipliers from Jacobian of `p`-step return map",
            "`scripts/dynamical_zeta/find_cycles.py` (`estimate_cycle_multiplier`)",
            (
                f"`unique_cycles={cycle['cycles_unique']}, first_period={cycle['first_period']}, "
                f"first_multiplier_abs={cycle['first_multiplier_abs']}, first_stability={cycle['first_stability']}` "
                f"(table: `{cycle['out_file']}`)"
            ),
        ),
        (
            "Gain/loss channel in driven systems",
            "anti-term yields real gain/loss factor after exponentiation",
            "`src/collatzscape/maps.py` (`phi_anti`, `c_eff`)",
            (
                f"`log|ce(-)|-log|ce(+)|={stroboscopic['log_ratio']:.12g}, "
                f"gamma*delta_f={stroboscopic['gamma_delta_f']:.12g}, "
                f"abs_error={stroboscopic['identity_abs_error']:.3e}`"
            ),
        ),
        (
            "Effective quasiparticle labels",
            "each `c` behaves as a parameterized mode with distinct survival/escape profile",
            "`figures/survival_map*.png`, sweep outputs",
            (
                "`"
                + ", ".join([f"c={k}:S={v:.4f}" for k, v in modes.items()])
                + "`"
            ),
        ),
    ]

    lines = [
        "# Floquet Connection Check",
        "",
        "| Floquet concept | Repo analog | Implemented evidence | Run result |",
        "|---|---|---|---|",
    ]
    for row in rows:
        lines.append(f"| {row[0]} | {row[1]} | {row[2]} | {row[3]} |")
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser(description="Run reproducible checks for Floquet-connection evidence table.")
    ap.add_argument("--c", default="1.8725+0j")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--gamma", type=float, default=0.15)
    ap.add_argument("--direction", type=int, choices=[-1, 1], default=1)
    ap.add_argument("--kick-mode", choices=["tanh", "clip"], default="tanh")
    ap.add_argument("--C-mode", choices=["cosine", "linear"], default="cosine")
    ap.add_argument("--z0", default="0.8+0.1j")
    ap.add_argument("--w0", default="0.85-0.1j")
    ap.add_argument("--out-report", default="outputs/csv/floquet_connection_check.md")
    ap.add_argument("--out-cycles", default="outputs/csv/dynamical_zeta/cycles_floquet_check.txt")
    args = ap.parse_args()

    c = complex(args.c.replace("i", "j"))
    z0 = complex(args.z0.replace("i", "j"))
    w0 = complex(args.w0.replace("i", "j"))
    out_report = Path(args.out_report).resolve()
    out_cycles = Path(args.out_cycles).resolve()

    stroboscopic = run_stroboscopic_and_gain_loss_check(
        c=c,
        alpha=args.alpha,
        gamma=args.gamma,
        direction=args.direction,
        z=z0,
        w=w0,
        kick_mode=args.kick_mode,
        C_mode=args.C_mode,
    )
    cycle = run_cycle_multiplier_check(
        c=c,
        alpha=args.alpha,
        gamma=args.gamma,
        direction=args.direction,
        kick_mode=args.kick_mode,
        C_mode=args.C_mode,
        out_path=out_cycles,
    )
    modes = run_quasiparticle_mode_check(alpha=args.alpha)

    report = build_markdown_report(stroboscopic=stroboscopic, cycle=cycle, modes=modes)
    out_report.parent.mkdir(parents=True, exist_ok=True)
    out_report.write_text(report + "\n", encoding="utf-8")

    print(report)
    print(f"\nSaved report: {out_report}")
    print(f"Saved cycle table: {out_cycles}")


if __name__ == "__main__":
    main()

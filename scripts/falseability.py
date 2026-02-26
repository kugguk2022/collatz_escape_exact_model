from __future__ import annotations

import argparse

import numpy as np


def C_cos_collatz(z: complex) -> complex:
    cz = np.cos(np.pi * z)
    return (z / 4.0) * (1.0 + cz) + ((3.0 * z + 1.0) / 2.0) * (1.0 - cz)


def phase_kick(cz: complex, alpha: float, phase_mode: str = "exp", phi_clip: float = 60.0, tanh_scale: float = 1.0) -> float:
    x = float(np.real(cz))
    if phase_mode == "exp":
        return alpha * float(np.exp(np.clip(x, -phi_clip, phi_clip)))
    if phase_mode == "tanh":
        return alpha * float(np.tanh(tanh_scale * x))
    raise ValueError(f"Unknown phase_mode: {phase_mode}")


def orbit_mean_g(
    c: complex,
    z0: complex,
    alpha: float,
    n_steps: int = 160,
    radius: float = 80.0,
    phase_mode: str = "exp",
    phi_clip: float = 60.0,
    tanh_scale: float = 1.0,
) -> tuple[bool, float, int]:
    """
    Returns (survived, mean_g, steps_used).
    survived=True if trajectory does not escape by n_steps.
    """
    z = np.complex128(z0)

    arg_c = np.angle(c)
    ln_abs_c = np.log(np.abs(c)) if np.abs(c) > 0 else -np.inf
    log_radius = np.log(radius)

    g_sum = 0.0
    steps = 0

    for _ in range(n_steps):
        cz = C_cos_collatz(z)

        # G_n proxy used in the original script.
        phi = phase_kick(cz, alpha=alpha, phase_mode=phase_mode, phi_clip=phi_clip, tanh_scale=tanh_scale)
        g_n = (np.imag(cz) * (arg_c + phi) - np.real(cz) * ln_abs_c)

        g_sum += float(np.real(g_n))
        steps += 1

        c_eff = c * np.exp(1j * phi)
        power_arg = cz * np.log(c_eff)
        if (not np.isfinite(power_arg.real)) or (not np.isfinite(power_arg.imag)):
            return (False, g_sum / steps, steps)
        if power_arg.real > log_radius:
            return (False, g_sum / steps, steps)

        z = np.exp(power_arg)

        if (not np.isfinite(z.real)) or (not np.isfinite(z.imag)) or (np.abs(z) > radius):
            return (False, g_sum / steps, steps)

    return (True, g_sum / steps, steps)


def compare_groups(
    c: complex,
    z0_mode: str = "c",
    alpha: float = 0.05,
    n_steps: int = 160,
    radius: float = 80.0,
    trials: int = 200,
    dz: float = 0.03,
    da: float = 0.01,
    seed: int = 0,
    phase_mode: str = "exp",
    phi_clip: float = 60.0,
    tanh_scale: float = 1.0,
) -> tuple[dict, dict]:
    rng = np.random.default_rng(seed)

    surv_g = []
    esc_g = []

    for _ in range(trials):
        z0_base = c if z0_mode == "c" else (1.0 + 0j)
        z0 = z0_base + (rng.normal(0, dz) + 1j * rng.normal(0, dz))
        a = alpha + rng.normal(0, da)

        survived, mean_g, _steps = orbit_mean_g(
            c,
            z0,
            a,
            n_steps=n_steps,
            radius=radius,
            phase_mode=phase_mode,
            phi_clip=phi_clip,
            tanh_scale=tanh_scale,
        )
        if survived:
            surv_g.append(mean_g)
        else:
            esc_g.append(mean_g)

    def summarize(xs: list[float]) -> dict:
        arr = np.array(xs, dtype=float)
        if arr.size == 0:
            return {"n": 0}
        arr = arr[np.isfinite(arr)]
        if arr.size == 0:
            return {"n": 0}
        return {
            "n": int(arr.size),
            "mean": float(arr.mean()),
            "median": float(np.median(arr)),
            "std": float(arr.std(ddof=1)) if arr.size > 1 else 0.0,
        }

    return summarize(surv_g), summarize(esc_g)


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare G-statistics for survivors vs escapers.")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--steps", type=int, default=160)
    parser.add_argument("--radius", type=float, default=80.0)
    parser.add_argument("--trials", type=int, default=300)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--dz", type=float, default=0.03, help="Jitter std for z0 real/imag.")
    parser.add_argument("--da", type=float, default=0.01, help="Jitter std for alpha.")
    parser.add_argument("--phase-mode", choices=["exp", "tanh"], default="exp")
    parser.add_argument("--tanh-scale", type=float, default=1.0)
    parser.add_argument("--deterministic", action="store_true", help="Disable jitter/noise (dz=da=0).")
    parser.add_argument(
        "--deterministic-tanh",
        action="store_true",
        help="Shortcut for deterministic tanh phase kick (phase-mode=tanh, dz=da=0).",
    )
    parser.add_argument("--phi-clip", type=float, default=60.0)
    parser.add_argument("--z0-mode", choices=["c", "one"], default="c")
    args = parser.parse_args()

    phase_mode = args.phase_mode
    dz = args.dz
    da = args.da
    if args.deterministic:
        dz = 0.0
        da = 0.0
    if args.deterministic_tanh:
        phase_mode = "tanh"
        dz = 0.0
        da = 0.0

    for c in [1.8725 + 0j, 2.20 + 0j, 1.87 + 0.35j]:
        surv_stats, esc_stats = compare_groups(
            c,
            z0_mode=args.z0_mode,
            alpha=args.alpha,
            n_steps=args.steps,
            radius=args.radius,
            trials=args.trials,
            dz=dz,
            da=da,
            seed=args.seed,
            phase_mode=phase_mode,
            phi_clip=args.phi_clip,
            tanh_scale=args.tanh_scale,
        )
        print(f"\nc = {c}")
        print(" survivors:", surv_stats)
        print(" escapers :", esc_stats)


if __name__ == "__main__":
    main()

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def C_cos_collatz(z: complex) -> complex:
    """Continuous cosine-interpolated Collatz gate."""
    cz = np.cos(np.pi * z)
    return (z / 4.0) * (1.0 + cz) + ((3.0 * z + 1.0) / 2.0) * (1.0 - cz)


def phase_kick(cz: complex, alpha: float, phase_mode: str = "exp", phi_clip: float = 60.0, tanh_scale: float = 1.0) -> float:
    """Phase kick variants used by this exploratory script."""
    x = float(np.real(cz))
    if phase_mode == "exp":
        return alpha * float(np.exp(np.clip(x, -phi_clip, phi_clip)))
    if phase_mode == "tanh":
        return alpha * float(np.tanh(tanh_scale * x))
    raise ValueError(f"Unknown phase_mode: {phase_mode}")


def survive_one(
    c: complex,
    z0: complex,
    alpha: float,
    n_steps: int,
    radius: float,
    phase_mode: str = "exp",
    phi_clip: float = 60.0,
    tanh_scale: float = 1.0,
) -> bool:
    """
    Returns True if trajectory survives to n_steps (never escapes and stays finite).
    Uses principal log branch.
    """
    z = np.complex128(z0)
    log_radius = np.log(radius)

    for _ in range(n_steps):
        cz = C_cos_collatz(z)

        phi = phase_kick(cz, alpha=alpha, phase_mode=phase_mode, phi_clip=phi_clip, tanh_scale=tanh_scale)
        c_eff = c * np.exp(1j * phi)

        power_arg = cz * np.log(c_eff)
        if (not np.isfinite(power_arg.real)) or (not np.isfinite(power_arg.imag)):
            return False

        # |exp(a+ib)| = exp(a); if Re(arg) > log(radius), fail fast.
        if power_arg.real > log_radius:
            return False

        z = np.exp(power_arg)

        if (not np.isfinite(z.real)) or (not np.isfinite(z.imag)) or (np.abs(z) > radius):
            return False

    return True


def survival_probability(
    c: complex,
    z0: complex,
    alpha: float,
    n_steps: int = 160,
    radius: float = 80.0,
    trials: int = 24,
    dz: float = 0.03,
    da: float = 0.01,
    seed: int = 0,
    phase_mode: str = "exp",
    phi_clip: float = 60.0,
    tanh_scale: float = 1.0,
) -> float:
    """
    S(c): fraction of jittered trials that survive.
    Jitter helps reduce one-lucky-orbit artifacts.
    """
    rng = np.random.default_rng(seed)
    ok = 0
    for _ in range(trials):
        z0_j = z0 + (rng.normal(0, dz) + 1j * rng.normal(0, dz))
        a_j = alpha + rng.normal(0, da)
        if survive_one(
            c,
            z0_j,
            a_j,
            n_steps,
            radius,
            phase_mode=phase_mode,
            phi_clip=phi_clip,
            tanh_scale=tanh_scale,
        ):
            ok += 1
    return ok / trials


def survival_map(
    re_min: float = -2.5,
    re_max: float = 1.5,
    im_min: float = -2.0,
    im_max: float = 2.0,
    resolution: int = 250,
    z0_mode: str = "c",
    alpha: float = 0.05,
    n_steps: int = 160,
    radius: float = 80.0,
    trials: int = 16,
    dz: float = 0.03,
    da: float = 0.01,
    seed: int = 0,
    phase_mode: str = "exp",
    phi_clip: float = 60.0,
    tanh_scale: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = np.linspace(re_min, re_max, resolution)
    y = np.linspace(im_min, im_max, resolution)
    grid_x, grid_y = np.meshgrid(x, y)
    c_grid = grid_x + 1j * grid_y

    surv = np.zeros(c_grid.shape, dtype=float)

    for iy in range(resolution):
        for ix in range(resolution):
            c = c_grid[iy, ix]
            if c == 0:
                surv[iy, ix] = 0.0
                continue
            z0 = c if z0_mode == "c" else (1.0 + 0j)
            surv[iy, ix] = survival_probability(
                c,
                z0=z0,
                alpha=alpha,
                n_steps=n_steps,
                radius=radius,
                trials=trials,
                dz=dz,
                da=da,
                seed=seed,
                phase_mode=phase_mode,
                phi_clip=phi_clip,
                tanh_scale=tanh_scale,
            )

    return grid_x, grid_y, surv


def render_map(grid_x: np.ndarray, grid_y: np.ndarray, surv: np.ndarray, output_png: str) -> None:
    fig = plt.figure(figsize=(8, 6))
    plt.pcolormesh(grid_x, grid_y, surv, shading="auto")
    plt.xlabel("Re(c)")
    plt.ylabel("Im(c)")
    plt.title("Survival probability S(c) for feedback-rotated Collatz-tetration")
    plt.colorbar(label="S(c) = fraction surviving to N")
    plt.grid(alpha=0.2)

    backend = plt.get_backend().lower()
    if "agg" in backend:
        out = Path(output_png).resolve()
        fig.savefig(out, dpi=180, bbox_inches="tight")
        print(f"Non-interactive backend '{plt.get_backend()}': saved {out}")
        plt.close(fig)
    else:
        plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(description="Quasiparticle-style survival map experiment.")
    parser.add_argument("--resolution", type=int, default=220)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--steps", type=int, default=160)
    parser.add_argument("--radius", type=float, default=80.0)
    parser.add_argument("--trials", type=int, default=12)
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
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--z0-mode", choices=["c", "one"], default="c")
    parser.add_argument("--output", default="figures/survival_map_quasiparticles.png")
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

    grid_x, grid_y, surv = survival_map(
        resolution=args.resolution,
        z0_mode=args.z0_mode,
        alpha=args.alpha,
        n_steps=args.steps,
        radius=args.radius,
        trials=args.trials,
        dz=dz,
        da=da,
        phi_clip=args.phi_clip,
        seed=args.seed,
        phase_mode=phase_mode,
        tanh_scale=args.tanh_scale,
    )
    render_map(grid_x, grid_y, surv, output_png=args.output)


if __name__ == "__main__":
    main()

from __future__ import annotations
import numpy as np
from .simulate import iterate
from .maps import step_pair

def escape_rate(points: np.ndarray, *, c: complex, alpha: float, gamma: float, direction: int, steps: int,
                kick_mode: str = "tanh", C_mode: str = "cosine", paired: bool = True) -> float:
    escaped = 0
    for z0 in points:
        _, esc, _ = iterate(complex(z0), c=c, alpha=alpha, gamma=gamma, direction=direction, steps=steps,
                            kick_mode=kick_mode, C_mode=C_mode, paired=paired)
        escaped += int(esc)
    return escaped / len(points)

def _step_unpaired_z(z: complex, *, c: complex, alpha: float, gamma: float, direction: int,
                     kick_mode: str, C_mode: str) -> complex:
    z_next, _ = step_pair(
        z,
        np.conjugate(z),
        c=c,
        alpha=alpha,
        gamma=gamma,
        direction=direction,
        kick_mode=kick_mode,
        C_mode=C_mode,
    )
    return z_next


def finite_time_lyapunov_rough(z0: complex, *, c: complex, alpha: float, gamma: float, direction: int,
                               steps: int, kick_mode: str = "tanh", C_mode: str = "cosine",
                               paired: bool = True, eps: float = 1e-7, warmup_steps: int = 0,
                               dt: float = 1.0, seed: int | None = 0) -> float:
    """
    Rough MLE estimate:
    - single Benettin trial (no ensemble averaging),
    - sensitive to perturbation direction and transients.
    """
    if steps <= 0:
        return float("nan")
    if eps <= 0.0:
        raise ValueError("eps must be > 0")
    if dt <= 0.0:
        raise ValueError("dt must be > 0")
    rng = np.random.default_rng(seed)
    return _mle_single_trial(
        z0,
        c=c,
        alpha=alpha,
        gamma=gamma,
        direction=direction,
        steps=steps,
        kick_mode=kick_mode,
        C_mode=C_mode,
        paired=paired,
        eps=eps,
        warmup_steps=warmup_steps,
        dt=dt,
        rng=rng,
    )


def _mle_single_trial(z0: complex, *, c: complex, alpha: float, gamma: float, direction: int, steps: int,
                      kick_mode: str, C_mode: str, paired: bool, eps: float, warmup_steps: int,
                      dt: float, rng: np.random.Generator) -> float:
    z_ref = complex(z0)
    w_ref = complex(np.conjugate(z0))

    for _ in range(max(0, int(warmup_steps))):
        if paired:
            z_ref, w_ref = step_pair(
                z_ref,
                w_ref,
                c=c,
                alpha=alpha,
                gamma=gamma,
                direction=direction,
                kick_mode=kick_mode,
                C_mode=C_mode,
            )
        else:
            z_ref = _step_unpaired_z(
                z_ref,
                c=c,
                alpha=alpha,
                gamma=gamma,
                direction=direction,
                kick_mode=kick_mode,
                C_mode=C_mode,
            )
            w_ref = np.conjugate(z_ref)

    if paired:
        delta = rng.normal(size=4)
        delta_norm = float(np.linalg.norm(delta))
        if delta_norm == 0.0:
            return float("nan")
        delta = (eps / delta_norm) * delta
        z_pert = z_ref + complex(delta[0], delta[1])
        w_pert = w_ref + complex(delta[2], delta[3])
    else:
        delta = rng.normal(size=2)
        delta_norm = float(np.linalg.norm(delta))
        if delta_norm == 0.0:
            return float("nan")
        delta = (eps / delta_norm) * delta
        z_pert = z_ref + complex(delta[0], delta[1])
        w_pert = np.conjugate(z_pert)

    log_sum = 0.0
    used_steps = 0

    for _ in range(steps):
        if paired:
            z_ref, w_ref = step_pair(
                z_ref,
                w_ref,
                c=c,
                alpha=alpha,
                gamma=gamma,
                direction=direction,
                kick_mode=kick_mode,
                C_mode=C_mode,
            )
            z_pert, w_pert = step_pair(
                z_pert,
                w_pert,
                c=c,
                alpha=alpha,
                gamma=gamma,
                direction=direction,
                kick_mode=kick_mode,
                C_mode=C_mode,
            )
            delta_vec = np.array(
                [z_pert.real - z_ref.real, z_pert.imag - z_ref.imag, w_pert.real - w_ref.real, w_pert.imag - w_ref.imag],
                dtype=float,
            )
        else:
            z_ref = _step_unpaired_z(
                z_ref,
                c=c,
                alpha=alpha,
                gamma=gamma,
                direction=direction,
                kick_mode=kick_mode,
                C_mode=C_mode,
            )
            w_ref = np.conjugate(z_ref)
            z_pert = _step_unpaired_z(
                z_pert,
                c=c,
                alpha=alpha,
                gamma=gamma,
                direction=direction,
                kick_mode=kick_mode,
                C_mode=C_mode,
            )
            w_pert = np.conjugate(z_pert)
            delta_vec = np.array([z_pert.real - z_ref.real, z_pert.imag - z_ref.imag], dtype=float)

        dist = float(np.linalg.norm(delta_vec))
        if (not np.isfinite(dist)) or dist <= 0.0:
            break

        log_sum += float(np.log(dist / eps))
        used_steps += 1

        scale = eps / dist
        z_pert = z_ref + (z_pert - z_ref) * scale
        if paired:
            w_pert = w_ref + (w_pert - w_ref) * scale
        else:
            w_pert = np.conjugate(z_pert)

    if used_steps == 0:
        return float("nan")
    return float(log_sum / (used_steps * dt))


def maximal_lyapunov_stats(z0: complex, *, c: complex, alpha: float, gamma: float, direction: int,
                           steps: int, kick_mode: str = "tanh", C_mode: str = "cosine", paired: bool = True,
                           eps: float = 1e-8, warmup_steps: int = 32, n_trials: int = 8,
                           seed: int | None = 0, dt: float = 1.0) -> dict[str, float | int]:
    """
    Benettin/Wolf-style maximal Lyapunov estimator on the map state.

    Returns:
        {"mle": mean, "mle_std": sample_std, "valid_trials": n}
    """
    if steps <= 0:
        return {"mle": float("nan"), "mle_std": float("nan"), "valid_trials": 0}
    if n_trials <= 0:
        raise ValueError("n_trials must be >= 1")
    if eps <= 0.0:
        raise ValueError("eps must be > 0")
    if dt <= 0.0:
        raise ValueError("dt must be > 0")

    rng = np.random.default_rng(seed)
    values: list[float] = []

    for _ in range(n_trials):
        val = _mle_single_trial(
            z0,
            c=c,
            alpha=alpha,
            gamma=gamma,
            direction=direction,
            steps=steps,
            kick_mode=kick_mode,
            C_mode=C_mode,
            paired=paired,
            eps=eps,
            warmup_steps=warmup_steps,
            dt=dt,
            rng=rng,
        )
        if np.isfinite(val):
            values.append(float(val))

    if not values:
        return {"mle": float("nan"), "mle_std": float("nan"), "valid_trials": 0}

    arr = np.array(values, dtype=float)
    std = float(np.std(arr, ddof=1)) if arr.size > 1 else 0.0
    return {"mle": float(np.mean(arr)), "mle_std": std, "valid_trials": int(arr.size)}


def finite_time_lyapunov(z0: complex, *, c: complex, alpha: float, gamma: float, direction: int,
                         steps: int, kick_mode: str = "tanh", C_mode: str = "cosine",
                         paired: bool = True, eps: float = 1e-8, warmup_steps: int = 32,
                         n_trials: int = 8, seed: int | None = 0, dt: float = 1.0) -> float:
    """
    Accurate maximal Lyapunov estimate (mean across randomized Benettin trials).
    """
    stats = maximal_lyapunov_stats(
        z0,
        c=c,
        alpha=alpha,
        gamma=gamma,
        direction=direction,
        steps=steps,
        kick_mode=kick_mode,
        C_mode=C_mode,
        paired=paired,
        eps=eps,
        warmup_steps=warmup_steps,
        n_trials=n_trials,
        seed=seed,
        dt=dt,
    )
    return float(stats["mle"])

from __future__ import annotations
import numpy as np
from .simulate import iterate

def escape_rate(points: np.ndarray, *, c: complex, alpha: float, gamma: float, direction: int, steps: int,
                kick_mode: str = "tanh", C_mode: str = "cosine", paired: bool = True) -> float:
    escaped = 0
    for z0 in points:
        _, esc, _ = iterate(complex(z0), c=c, alpha=alpha, gamma=gamma, direction=direction, steps=steps,
                            kick_mode=kick_mode, C_mode=C_mode, paired=paired)
        escaped += int(esc)
    return escaped / len(points)

def finite_time_lyapunov(z0: complex, *, c: complex, alpha: float, gamma: float, direction: int,
                         steps: int, kick_mode: str = "tanh", C_mode: str = "cosine",
                         paired: bool = True, eps: float = 1e-7) -> float:
    """
    Rough finite-time Lyapunov estimate via two nearby trajectories in z.
    """
    z = z0
    z2 = z0 + eps
    for _ in range(steps):
        z_next = iterate(z, c=c, alpha=alpha, gamma=gamma, direction=direction, steps=1,
                         kick_mode=kick_mode, C_mode=C_mode, paired=paired)[0]
        z2_next = iterate(z2, c=c, alpha=alpha, gamma=gamma, direction=direction, steps=1,
                          kick_mode=kick_mode, C_mode=C_mode, paired=paired)[0]
        dz = abs(z2_next - z_next)
        if dz == 0 or (not np.isfinite(dz)):
            return float("nan")
        z2 = z_next + (z2_next - z_next) * (eps / dz)
        z = z_next
    return float(np.log(abs((z2 - z) / eps)))

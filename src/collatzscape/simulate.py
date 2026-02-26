from __future__ import annotations
import numpy as np
from .maps import step_pair

def iterate(z0: complex, *, c: complex, alpha: float, gamma: float, direction: int, steps: int,
            kick_mode: str = "tanh", C_mode: str = "cosine",
            escape_radius: float = 1e6, paired: bool = True) -> tuple[complex, bool, int]:
    """
    Iterate starting at z0.

    paired=True: evolve (z,w) with w0=conj(z0).
    paired=False: use w=conj(z) each step (no independent w evolution).
    """
    z = complex(z0)
    w = complex(np.conjugate(z0))

    for t in range(steps):
        if paired:
            z, w = step_pair(z, w, c=c, alpha=alpha, gamma=gamma, direction=direction,
                             kick_mode=kick_mode, C_mode=C_mode)
        else:
            z, _ = step_pair(z, np.conjugate(z), c=c, alpha=alpha, gamma=gamma, direction=direction,
                             kick_mode=kick_mode, C_mode=C_mode)
            w = np.conjugate(z)

        # Guard against huge magnitudes: non-finite values and overflow both count as escape.
        if (not np.isfinite(z.real)) or (not np.isfinite(z.imag)):
            return z, True, t
        try:
            if abs(z) > escape_radius:
                return z, True, t
        except OverflowError:
            return z, True, t
    return z, False, steps

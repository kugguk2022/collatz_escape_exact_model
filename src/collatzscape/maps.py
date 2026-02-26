"""
maps.py — Cosine-Collatz gate + phase-kick framework (repo-grounded)

Core update:
    z_{n+1} = exp( C(z_n) * Log( c_eff(z_n, w_n) ) )

Cosine interpolant (entire, complex-extendable):
    C(z) = (z/4)*(1 + cos(pi z)) + ((3z + 1)/2)*(1 - cos(pi z))

Phase kick (paired + anti):
    f(u) = tanh(u)   (or clip)
    phi_base(z) = alpha * f(Re(C(z)))
    phi_paired(z,w) = 0.5*(phi_base(z) + phi_base(w))
    phi_anti(z,w) = direction*(1j*gamma/2)*( f(Re(C(z))) - f(Re(C(w))) )

    phi_total = phi_paired + phi_anti
    c_eff = c * exp( 1j * phi_total )

Notes:
- If w == conj(z) exactly and C respects conjugation (it does), then Re(C(w)) == Re(C(z)),
  hence phi_anti = 0 and phi_total is real → |c_eff| = |c| (pure phase).
- If w only approximately tracks conj(z), phi_anti can become nonzero and can introduce
  direction-dependent damping/amplification via exp(1j*(a+ib)) = exp(i a - b).
"""
from __future__ import annotations
import numpy as np

def C_cosine(z: complex) -> complex:
    """Cosine-interpolated Collatz gate extended to complex z (entire)."""
    return (z/4.0)*(1.0 + np.cos(np.pi*z)) + ((3.0*z + 1.0)/2.0)*(1.0 - np.cos(np.pi*z))

def C_linear_near2(z: complex) -> complex:
    """Local linearization around z=2: C(z)=1 + (1/2)(z-2) + ..."""
    return 1.0 + 0.5*(z - 2.0)

def f_tanh(u: float) -> float:
    """Bounded in [-1,1]."""
    return float(np.tanh(u))

def f_clip(u: float, lo: float = -1.0, hi: float = 1.0) -> float:
    """Bounded in [lo,hi]."""
    return float(np.clip(u, lo, hi))

def phi_base_from_ReC(ReC: float, *, alpha: float, mode: str = "tanh") -> float:
    """
    Real-valued phase-kick amplitude.
    Since tanh/clip outputs in [-1,1], the effective kick range is [-alpha, +alpha].
    """
    if mode == "tanh":
        return alpha * f_tanh(ReC)
    if mode == "clip":
        return alpha * f_clip(ReC)
    raise ValueError(f"Unknown kick mode: {mode}")

def c_eff(z: complex, w: complex, *, c: complex, alpha: float, gamma: float,
          direction: int = +1, kick_mode: str = "tanh", C_mode: str = "cosine") -> complex:
    """Compute c_eff(z,w)."""
    Cfun = C_cosine if C_mode == "cosine" else C_linear_near2
    Cz = Cfun(z)
    Cw = Cfun(w)

    pb_z = phi_base_from_ReC(float(np.real(Cz)), alpha=alpha, mode=kick_mode)
    pb_w = phi_base_from_ReC(float(np.real(Cw)), alpha=alpha, mode=kick_mode)

    phi_paired = 0.5*(pb_z + pb_w)  # real

    sig_z = f_tanh(float(np.real(Cz))) if kick_mode == "tanh" else f_clip(float(np.real(Cz)))
    sig_w = f_tanh(float(np.real(Cw))) if kick_mode == "tanh" else f_clip(float(np.real(Cw)))
    phi_anti = direction * (1j * gamma / 2.0) * (sig_z - sig_w)

    phi_total = phi_paired + phi_anti
    return c * np.exp(1j * phi_total)

def step_pair(z: complex, w: complex, *, c: complex, alpha: float = 0.0, gamma: float = 0.0,
              direction: int = +1, kick_mode: str = "tanh", C_mode: str = "cosine") -> tuple[complex, complex]:
    """
    One full step updating the paired variables (z, w):

        z_{n+1} = exp( C(z_n) * Log( c_eff(z_n, w_n) ) )
        w_{n+1} = exp( C(w_n) * Log( c_eff(w_n, z_n) ) )
    """
    Cfun = C_cosine if C_mode == "cosine" else C_linear_near2
    Cz = Cfun(z)
    Cw = Cfun(w)

    ce_z = c_eff(z, w, c=c, alpha=alpha, gamma=gamma, direction=direction, kick_mode=kick_mode, C_mode=C_mode)
    ce_w = c_eff(w, z, c=c, alpha=alpha, gamma=gamma, direction=direction, kick_mode=kick_mode, C_mode=C_mode)

    z_next = np.exp(Cz * np.log(ce_z))
    w_next = np.exp(Cw * np.log(ce_w))
    return complex(z_next), complex(w_next)

def step(z: complex, *, c: complex, alpha: float = 0.0, gamma: float = 0.0, direction: int = +1,
         kick_mode: str = "tanh", C_mode: str = "cosine", w: complex | None = None) -> complex:
    """Convenience wrapper returning z_{n+1}. If w is None, uses w=conj(z) for that step."""
    if w is None:
        w = np.conjugate(z)
    z_next, _ = step_pair(z, w, c=c, alpha=alpha, gamma=gamma, direction=direction, kick_mode=kick_mode, C_mode=C_mode)
    return z_next

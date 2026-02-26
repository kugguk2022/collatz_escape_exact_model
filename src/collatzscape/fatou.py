from __future__ import annotations
import numpy as np
from .simulate import iterate
from .maps import step

def jacobian_real(z: complex, *, c: complex, gamma: float, direction: int, h: float = 1e-7) -> np.ndarray:
    """
    Numerical Jacobian of the 2D real map induced by the complex update.

    Treat z = x + i y, map T: R^2 -> R^2.
    """
    x, y = float(np.real(z)), float(np.imag(z))

    def T(xy):
        zz = complex(float(xy[0]), float(xy[1]))
        out = step(zz, c=c, gamma=gamma, direction=direction)
        return np.array([np.real(out), np.imag(out)], dtype=float)

    p = np.array([x, y], dtype=float)
    e0 = np.array([h, 0.0], dtype=float)
    e1 = np.array([0.0, h], dtype=float)
    f0 = T(p)
    f1 = T(p + e0)
    f2 = T(p + e1)
    J = np.column_stack(((f1 - f0) / h, (f2 - f0) / h))
    return J

def local_stability(z_star: complex, *, c: complex, gamma: float, direction: int) -> dict:
    """
    Return eigenvalues and spectral radius of the local linearization around z_star.
    For attracting behavior in the real 2D sense, we want spectral_radius < 1.
    """
    J = jacobian_real(z_star, c=c, gamma=gamma, direction=direction)
    eig = np.linalg.eigvals(J)
    rho = float(np.max(np.abs(eig)))
    return {"eig1": complex(eig[0]), "eig2": complex(eig[1]), "spectral_radius": rho}

def find_attractors(*, c: complex, gamma: float, direction: int, steps: int,
                    re_min: float, re_max: float, im_min: float, im_max: float,
                    n_samples: int = 5000, tol: float = 1e-3, seed: int = 0) -> list[complex]:
    """
    Heuristic attractor finder:
    - sample random initial conditions in a window
    - iterate
    - cluster non-escaping endpoints by rounding to 'tol'
    Returns a list of representative points.
    """
    rng = np.random.default_rng(seed)
    xs = rng.uniform(re_min, re_max, size=n_samples)
    ys = rng.uniform(im_min, im_max, size=n_samples)
    reps = {}
    for x, y in zip(xs, ys):
        z0 = complex(float(x), float(y))
        zf, esc, _ = iterate(z0, c=c, gamma=gamma, direction=direction, steps=steps)
        if esc:
            continue
        key = (round(np.real(zf)/tol)*tol, round(np.imag(zf)/tol)*tol)
        reps.setdefault(key, 0)
        reps[key] += 1

    # sort by support count
    items = sorted(reps.items(), key=lambda kv: kv[1], reverse=True)
    return [complex(k[0], k[1]) for k, _ in items]

def classify_point(z0: complex, attractors: list[complex], *, c: complex, gamma: float, direction: int,
                   steps: int, escape_radius: float = 1e6) -> int:
    """
    Returns:
      -1 if escaped
      i (0..len(attractors)-1) for nearest attractor representative by endpoint.
    """
    zf, esc, _ = iterate(z0, c=c, gamma=gamma, direction=direction, steps=steps, escape_radius=escape_radius)
    if esc or len(attractors) == 0:
        return -1
    d = [abs(zf - a) for a in attractors]
    return int(np.argmin(d))

def fatou_label_grid(*, re_min: float, re_max: float, im_min: float, im_max: float,
                     n_re: int, n_im: int, attractors: list[complex],
                     c: complex, gamma: float, direction: int, steps: int,
                     escape_radius: float = 1e6) -> np.ndarray:
    """
    Grid labels:
      -1 escape
      0..K-1 basin labels (proxy for Fatou components / basins of attraction)
    """
    xs = np.linspace(re_min, re_max, n_re)
    ys = np.linspace(im_min, im_max, n_im)
    lab = np.full((n_im, n_re), -1, dtype=int)
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            z0 = complex(float(x), float(y))
            lab[j, i] = classify_point(z0, attractors, c=c, gamma=gamma, direction=direction,
                                       steps=steps, escape_radius=escape_radius)
    return lab

def basin_measure(labels: np.ndarray) -> dict:
    """
    Monte-carlo/grid proxy for basin measures.
    """
    total = labels.size
    escaped = int(np.sum(labels < 0))
    out = {"escaped": escaped/total}
    ks = sorted(set(int(k) for k in np.unique(labels) if k >= 0))
    for k in ks:
        out[f"basin_{k}"] = int(np.sum(labels == k))/total
    return out

"""
IEEE 14-bus: Collatz–Kuramoto gate + "power-tower" cascade risk (toy dynamics).

- Network data: IEEE 14-bus MATPOWER case14 (bus/gen/branch embedded).
- Dynamics: Kuramoto-with-inertia (swing-like) on line couplings K_ij ~ 1/(x*tap).
- Collatz gate: per-bus integer severity n_i updated from stress s_i then one Collatz step
    even -> n/2, odd -> 3n+1
  odd parity => "stressed" mode: more damping + weaker coupling
- Protection-like trip: if overloaded AND endpoint stressed => trip line (remove coupling)
- "Power tower" risk: stabilized nested amplification of (frequency, overload, incoherence)

This is a *concept demo*; it’s meant to be modified toward real protection/control logic.
"""
from __future__ import annotations
import numpy as np
import math

BASE_MVA = 100.0
BUS = np.array([[1, 3, 0, 0, 0, 0, 1, 1.06, 0, 0, 1, 1.06, 0.94], [2, 2, 21.7, 12.7, 0, 0, 1, 1.045, -4.98, 0, 1, 1.06, 0.94], [3, 2, 94.2, 19, 0, 0, 1, 1.01, -12.72, 0, 1, 1.06, 0.94], [4, 1, 47.8, -3.9, 0, 0, 1, 1.019, -10.33, 0, 1, 1.06, 0.94], [5, 1, 7.6, 1.6, 0, 0, 1, 1.02, -8.78, 0, 1, 1.06, 0.94], [6, 2, 11.2, 7.5, 0, 0, 1, 1.07, -14.22, 0, 1, 1.06, 0.94], [7, 1, 0, 0, 0, 0, 1, 1.062, -13.37, 0, 1, 1.06, 0.94], [8, 2, 0, 0, 0, 0, 1, 1.09, -13.36, 0, 1, 1.06, 0.94], [9, 1, 29.5, 16.6, 0, 19, 1, 1.056, -14.94, 0, 1, 1.06, 0.94], [10, 1, 9, 5.8, 0, 0, 1, 1.051, -15.1, 0, 1, 1.06, 0.94], [11, 1, 3.5, 1.8, 0, 0, 1, 1.057, -14.79, 0, 1, 1.06, 0.94], [12, 1, 6.1, 1.6, 0, 0, 1, 1.055, -15.07, 0, 1, 1.06, 0.94], [13, 1, 13.5, 5.8, 0, 0, 1, 1.05, -15.16, 0, 1, 1.06, 0.94], [14, 1, 14.9, 5, 0, 0, 1, 1.036, -16.04, 0, 1, 1.06, 0.94]], dtype=float)
GEN = np.array([[1, 232.4, -16.9, 10, 0, 1.06, 100, 1, 332.4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [2, 40, 42.4, 50, -40, 1.045, 100, 1, 140, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [3, 0, 23.4, 40, 0, 1.01, 100, 1, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [6, 0, 12.2, 24, -6, 1.07, 100, 1, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [8, 0, 17.4, 24, -6, 1.09, 100, 1, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype=float)
BR  = np.array([[1, 2, 0.01938, 0.05917, 0.0528, 1, 0, 0, 0, 0, 1, -60, 60], [1, 5, 0.05403, 0.22304, 0.0492, 489, 0, 0, 0, 0, 1, -60, 60], [2, 3, 0.04699, 0.19797, 0.0438, 552, 0, 0, 0, 0, 1, -60, 60], [2, 4, 0.05811, 0.17632, 0.034, 605, 0, 0, 0, 0, 1, -60, 60], [2, 5, 0.05695, 0.17388, 0.0346, 614, 0, 0, 0, 0, 1, -60, 60], [3, 4, 0.06701, 0.17103, 0.0128, 611, 0, 0, 0, 0, 1, -60, 60], [4, 5, 0.01335, 0.04211, 0, 2543, 0, 0, 0, 0, 1, -60, 60], [4, 7, 0, 0.20912, 0, 537, 0, 0, 0.978, 0, 1, -60, 60], [4, 9, 0, 0.55618, 0, 202, 0, 0, 0.969, 0, 1, -60, 60], [5, 6, 0, 0.25202, 0, 445, 0, 0, 0.932, 0, 1, -60, 60], [6, 11, 0.09498, 0.1989, 0, 509, 0, 0, 0, 0, 1, -60, 60], [6, 12, 0.12291, 0.25581, 0, 395, 0, 0, 0, 0, 1, -60, 60], [6, 13, 0.06615, 0.13027, 0, 769, 0, 0, 0, 0, 1, -60, 60], [7, 8, 0, 0.17615, 0, 637, 0, 0, 0, 0, 1, -60, 60], [7, 9, 0, 0.11001, 0, 1021, 0, 0, 0, 0, 1, -60, 60], [9, 10, 0.03181, 0.0845, 0, 1244, 0, 0, 0, 0, 1, -60, 60], [9, 14, 0.12711, 0.27038, 0, 376, 0, 0, 0, 0, 1, -60, 60], [10, 11, 0.08205, 0.19207, 0, 537, 0, 0, 0, 0, 1, -60, 60], [12, 13, 0.22092, 0.19988, 0, 377, 0, 0, 0, 0, 1, -60, 60], [13, 14, 0.17093, 0.34802, 0, 289, 0, 0, 0, 0, 1, -60, 60]], dtype=float)

def build_coupling(branch: np.ndarray):
    nb = int(np.max(branch[:, :2]))
    K = np.zeros((nb, nb), dtype=float)
    rateA = np.zeros((nb, nb), dtype=float)
    for row in branch:
        f, t = int(row[0])-1, int(row[1])-1
        x = float(row[3])
        tap = float(row[8]) if float(row[8]) != 0.0 else 1.0
        if x == 0.0:
            continue
        kij = 1.0/(x*tap)          # DC-ish coupling proxy
        K[f, t] += kij
        K[t, f] += kij
        ra = float(row[5])
        if ra <= 1.0:
            ra = 250.0
        rateA[f, t] = rateA[t, f] = ra
    return K, rateA

Kij, rateA = build_coupling(BR)

def simulate(
    T: float = 12.0,
    dt: float = 0.01,
    seed: int = 2,
    inertia_gen: float = 6.0,
    inertia_load: float = 1.5,
    damping_gen: float = 1.0,
    damping_load: float = 1.2,
    omega_th: float = 0.5,
    overload_trip: float = 1.35,
):
    rng = np.random.default_rng(seed)
    nb = BUS.shape[0]
    is_gen_bus = np.isin(np.arange(1, nb+1), GEN[:,0].astype(int))

    # states: angles (rad), frequency deviation proxy ω
    theta = np.deg2rad(BUS[:,8])
    omega = np.zeros(nb, dtype=float)

    # net injections Pm = Pg - Pd (per-unit)
    Pd = BUS[:,2] / BASE_MVA
    Pg = np.zeros(nb, dtype=float)
    for g in GEN:
        Pg[int(g[0])-1] += g[1] / BASE_MVA
    Pm = Pg - Pd + 0.02*rng.normal(size=nb)

    M = np.where(is_gen_bus, inertia_gen, inertia_load)
    D = np.where(is_gen_bus, damping_gen, damping_load)

    # line status (symmetric)
    status = (Kij > 0)

    steps = int(T/dt)
    out = {
        "t": np.linspace(0, T, steps),
        "r": np.zeros(steps),
        "max_over": np.zeros(steps),
        "max_omega": np.zeros(steps),
        "risk": np.zeros(steps),
        "tripped_edges": np.zeros(steps, dtype=int),
    }

    n = np.ones(nb, dtype=np.int64)

    for k in range(steps):
        out["r"][k] = np.abs(np.mean(np.exp(1j*theta)))

        dth = theta[:,None] - theta[None,:]
        Pline = np.zeros((nb, nb), dtype=float)
        Pline[status] = Kij[status] * dth[status]
        over = np.zeros((nb, nb), dtype=float)
        over[status] = np.abs(Pline[status]) / rateA[status]
        out["max_over"][k] = float(np.max(over)) if np.any(status) else 0.0

        inc_over = np.sum(over, axis=1)
        s = np.abs(omega)/omega_th + 0.5*inc_over

        # Collatz gate (integer severity)
        n = np.maximum(1, (1 + np.floor(10*s)).astype(np.int64))
        even = (n % 2 == 0)
        n_next = n.copy()
        n_next[even] = n[even]//2
        n_next[~even] = np.minimum(3*n[~even] + 1, 10_000_000)
        n = n_next
        stressed = (~even)

        # stressed mode: more damping + weaker coupling around stressed nodes
        D_eff = D * (1.0 + 1.0*stressed.astype(float))
        K_eff = Kij.copy()
        for i in range(nb):
            if stressed[i]:
                K_eff[i,:] *= 0.55
                K_eff[:,i] *= 0.55

        # protection-like trip
        trip = status & (over > overload_trip) & (stressed[:,None] | stressed[None,:])
        if np.any(trip):
            status[trip] = False
            # enforce symmetry
            ii, jj = np.nonzero(trip)
            status[jj, ii] = False

        out["tripped_edges"][k] = int(np.sum(status == False) - np.sum(Kij == 0))

        # electrical power term (sin coupling)
        sdiff = np.sin(theta[:,None] - theta[None,:])
        Pe = np.sum((K_eff * status) * sdiff, axis=1)

        domega = (Pm - Pe - D_eff*omega) / M
        omega = omega + dt*domega
        theta = theta + dt*omega

        out["max_omega"][k] = float(np.max(np.abs(omega)))

        # tower-style risk (stabilized)
        s_freq = min(1.0, out["max_omega"][k]/1.0)
        s_line = min(1.0, out["max_over"][k]/2.0)
        s_coh  = min(1.0, 1.0 - out["r"][k])
        b1 = math.exp(2.5*s_freq); b2 = math.exp(2.0*s_line); b3 = math.exp(1.5*s_coh)
        tower = math.log(b1) * (b2**(0.6 + 0.4*b3))
        out["risk"][k] = math.tanh(tower/8.0)

    return out

if __name__ == "__main__":
    out = simulate()
    print("final:", {
        "risk": float(out["risk"][-1]),
        "max_over": float(out["max_over"][-1]),
        "max_omega": float(out["max_omega"][-1]),
        "r": float(out["r"][-1]),
    })

# Floquet-Particle Connection: Similarities, Differences, Novelty

This note connects the current repository to Floquet-quasiparticle language while separating:
- what is directly analogous,
- what is different,
- what is genuinely new here.

## 1. Conservative Position

Collatzscape is best described as a **Floquet-like nonlinear stroboscopic map**, not a standard linear periodically driven Floquet Hamiltonian.

## 2. What Matches Floquet Intuition

| Floquet concept | Repo analog | Implemented evidence |
|---|---|---|
| Stroboscopic update over one drive period | one-step map `(z_n,w_n) -> (z_{n+1},w_{n+1})` | `src/collatzscape/maps.py` (`step_pair`) |
| Period-`p` orbit stability via multipliers | cycle multipliers from Jacobian of `p`-step return map | `scripts/dynamical_zeta/find_cycles.py` (`estimate_cycle_multiplier`) |
| Gain/loss channel in driven systems | anti-term yields real gain/loss factor after exponentiation | `src/collatzscape/maps.py` (`phi_anti`, `c_eff`) |
| Effective quasiparticle labels | each `c` behaves as a parameterized mode with distinct survival/escape profile | `figures/survival_map*.png`, sweep outputs |

Key anti-term identity (implemented):

\[
\phi_{\text{anti}}=\mathrm{direction}\cdot\frac{i\gamma}{2}\Delta_f,\quad
\Delta_f=f(\Re(C(z)))-f(\Re(C(w))).
\]

Inside `c_eff = c * exp(i * phi_total)`, this contributes:

\[
\exp(i\phi_{\text{anti}})=\exp\!\left(-\mathrm{direction}\cdot\frac{\gamma}{2}\Delta_f\right),
\]

which is an explicit directional damping/amplification factor.

## 3. Important Differences from Canonical Floquet Particles

1. **No external periodic drive is required**:
the forcing is endogenous/state-dependent through `Re(C(z_n))`, not a fixed-time periodic Hamiltonian.

2. **Strongly nonlinear map, not linear unitary evolution**:
the core update is transcendental (`exp(C(z) * Log(c_eff))`), branch-sensitive, and can overflow/escape.

3. **No quasi-energy band structure is currently defined**:
the repo computes escape rates, Lyapunov exponents, multipliers, and zeta-like products, not Bloch/Floquet bands.

4. **“Particles” are computational quasiparticles**:
they are parameter-space dynamical modes, not experimentally validated condensed-matter quasiparticles.

## 4. Novelty in This Repo (Relative to Standard Floquet Framing)

1. **Direction-odd anti-coupling with exact conjugate-manifold cancellation**:
if `w = conj(z)` exactly, `Delta_f = 0` and anti-channel vanishes.

2. **Exact diode observable for online gamma estimation**:
`y = log|c_eff(dir=-1)| - log|c_eff(dir=+1)| = gamma * Delta_f`
implemented and used in closed-loop control (`scripts/braid_diod.py`).

3. **Type-2 sector invariants loop**:
the run couples:
- tower fixed-point solve `x = a^x`,
- prime-sector `kappa_t`,
- odd-harmonic (Catalan-weighted) energy `gamma_{t+1}`,
into a single adaptive dynamics (`scripts/Type_2_sector_invariants.py`).

4. **Single pipeline across map-level and cycle-level diagnostics**:
escape/Lyapunov/basin maps plus cycle multipliers and finite-product zeta scans coexist in one reproducible codebase.

## 5. Where Pisot Relations Stand

A Pisot link is currently a **hypothesis/future direction**, not an implemented validated result in this snapshot.
Recommended phrasing: “Pisot-related structure is an open conjectural connection for future testing.”

## 6. Suggested Paper Wording (Safe)

“We interpret Collatzscape as a Floquet-like nonlinear stroboscopic system: period-`p` cycle multipliers play the role of local Floquet multipliers, while an endogenous anti-coupling term induces direction-selective gain/loss. Unlike canonical Floquet Hamiltonians, the drive is state-dependent and the dynamics are branch-sensitive transcendental iterates rather than linear unitary evolution.”

# About Collatzscape

Collatzscape is a reproducible framework for a forced complex map with direction-dependent stability.
The model combines a Collatz-style continuous gate with exponential/tetration-like iteration.

## Explicit Grounding

The project is grounded in non-autonomous/skew-product complex dynamics:

- map iteration with endogenous forcing,
- basin versus escape behavior,
- finite-time stability diagnostics (escape rate and Lyapunov proxies).

Canonical equations are in `MODEL.md` and `docs/00_foundations.md`.

## Core Formula (Short Version)

\[
$z_{n+1}=\exp\!\left(C(z_n)\Log(c_{\text{eff}}(z_n,w_n))\right)$,
\]

with \(C(z)\) the cosine-interpolated Collatz gate and $\(c_{\text{eff}}\)$ including paired and anti coupling.

## Applications

1. Benchmarking direction-dependent transitions in forced transcendental maps.
2. Building reproducible basin/escape atlases over parameter sweeps.
3. Stress-testing numerical handling of branch-sensitive exponential dynamics.
4. Rapid prototyping for deterministic versus jittered phase-kick regimes.
5. Exploratory analogy layer (quasiparticle/fusion language) for interpretive modeling.

## Scope and Boundaries

- This is a computational dynamics model, not a claim of a new Collatz theorem.
- Quasiparticle/fusion interpretation is optional and analogical.
- Main claim: explicit equations + reproducible experiments + clear diagnostics.

## Exploratory Appendix: Fusion-Tree Analogy

An optional interpretation maps paired states into an Ising-anyon-style fusion tree:

- total fusion to vacuum channel (`1`) as a global consistency condition,
- intermediate `1`/`psi` channel as a logical two-state subspace,
- braid-like operations interpreted as channel transformations.

This appendix is included for exploratory context and is not required for the core dynamical results.

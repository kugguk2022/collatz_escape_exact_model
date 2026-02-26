# About Collatzscape

Collatzescape is a reproducible research repo for a non-autonomous complex-dynamics system:
a Collatz-like gate drives an exponential/tetration-style update, and direction-dependent coupling
produces asymmetric stability (bounded behavior in one direction, escape in the other).

This file captures the exploratory interpretation from the shared notes linked by the project author.

## Fusion-Tree Interpretation (Three-Particle T-Junction)

In the paired model, three quasiparticle-like modes `(c1, c2, c3)` are treated with an Ising-anyon style fusion tree.

1. **Bottom layer (total fusion)**  
   The full three-particle system must fuse to vacuum (`1`) to preserve total topological charge.
   In this model, approximate pairing symmetry (`z_i ≈ conj(w_i)`) plays that role; strong pairing break
   corresponds to escape/loss of the vacuum-fusion condition.

2. **Intermediate fusion channel (logical qubit)**  
   The pair `(c1, c2)` fuses into either:
   - `1` (even parity, vacuum channel) -> logical `|0>`
   - `psi` (odd parity, fermion channel) -> logical `|1>`
   This 2D channel space is the protected logical subspace.

3. **Top layer (braiding input)**  
   Braiding `c1` around `c2` acts on that intermediate channel and changes the encoded logical state
   without a direct local measurement.

## Effective Braid Action

In the `{1, psi}` basis, the forward braid can be written (up to global phase) as:

```text
R = diag(e^{i*pi/8}, e^{i*5*pi/8})
```

This is the non-Abelian gate-like part of the interpretation: braiding acts on channel state, not just position.

## Directional Collapse Observation

From the exploratory sweeps in the shared content:

- Forward direction retains much higher logical-channel stability.
- Reverse direction can show a sharp collapse near a critical anti-term coupling (`gamma ~ 0.5` in that run).
- At that critical point, intermediate-channel coherence degrades strongly (fusion outcome randomization / leakage).

## Scope

This is an exploratory dynamical analogy built on the Collatzscape map, not a claim of physical Majorana hardware.
Treat it as a model-layer interpretation and a guide for further simulation/testing.

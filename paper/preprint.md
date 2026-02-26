# Collatzscape (repo-first preprint)

> This is a repo-first “paper” in Markdown. Keep it lightweight; regenerate figures from code.

## Abstract
We present an open implementation of a forced complex map with Collatz-style gating and exponential update. The system exhibits direction-dependent stability: one direction can remain bounded while the opposite direction develops dominant escape. We provide explicit equations, reproducible sweeps, and diagnostics (escape rates, finite-time Lyapunov proxies, and basin-component proxies).

## 1. Setup

Gate:

\[
C(z)=\frac{z}{4}(1+\cos(\pi z))+\frac{3z+1}{2}(1-\cos(\pi z)).
\]

Coupling:

\[
\phi_{\text{base}}(z)=\alpha f(\Re(C(z))),\quad
\phi_{\text{paired}}=\frac12(\phi_{\text{base}}(z)+\phi_{\text{base}}(w)),
\]
\[
\phi_{\text{anti}}=\mathrm{direction}\cdot\frac{i\gamma}{2}\left[f(\Re(C(z)))-f(\Re(C(w)))\right],
\]
\[
c_{\text{eff}}=c\exp(i(\phi_{\text{paired}}+\phi_{\text{anti}})).
\]

Update:

\[
z_{n+1}=\exp(C(z_n)\Log(c_{\text{eff}}(z_n,w_n))),
\quad
w_{n+1}=\exp(C(w_n)\Log(c_{\text{eff}}(w_n,z_n))).
\]

## 2. Direction-Dependent Transition

Because:

\[
\exp(i\phi_{\text{anti}})
=\exp\!\left(-\mathrm{direction}\cdot\frac{\gamma}{2}\Delta_f\right),
\quad
\Delta_f=f(\Re(C(z)))-f(\Re(C(w))),
\]

flipping `direction` flips effective gain/loss under mismatch \(\Delta_f\neq0\). This creates asymmetric collapse/escape channels.

## 3. Diagnostics

- Escape-rate estimator over sampled seeds.
- Finite-time Lyapunov proxy from nearby trajectories.
- Basin-component proxy from grid labeling and attractor clustering.
- Optional deterministic tanh/no-noise script variants for stress tests.

## 4. Applications and Limits

Applications:

1. Benchmarking non-reciprocal stability transitions.
2. Basin atlas generation in forced transcendental maps.
3. Numerical robustness testing for branch-sensitive exponentiation.

Limits:

- Not a proof of a new Collatz theorem.
- Quasiparticle/fusion language is optional analogy, not physical-device validation.

## 5. Reproducibility

All core figures are generated from pinned configs in `configs/` via:

```bash
python -m collatzscape.cli demo --config configs/default.yaml --out figures/
python -m collatzscape.cli sweep --config configs/default.yaml --out figures/
python -m collatzscape.cli fatou --config configs/default.yaml --out figures/
```

## References

See `docs/references.md`.

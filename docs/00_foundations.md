# Foundations: Grounding, Model, and Applications

## 1. Scope

Collatzscape studies a forced complex map where a Collatz-style gate modulates exponential iteration.
The goal is not to claim a new theorem, but to provide a reproducible model where directional stability asymmetry can be measured.

## 2. Canonical Map (Implemented)

Cosine-interpolated Collatz gate:

\[
C(z)=\frac{z}{4}\bigl(1+\cos(\pi z)\bigr)+\frac{3z+1}{2}\bigl(1-\cos(\pi z)\bigr).
\]

Bounded signal (default):

\[
f(u)=\tanh(u), \quad f(u)\in[-1,1].
\]

Base and coupling terms:

\[
\phi_{\text{base}}(z)=\alpha f(\Re(C(z))),
\]
\[
\phi_{\text{paired}}(z,w)=\frac12\left(\phi_{\text{base}}(z)+\phi_{\text{base}}(w)\right),
\]
\[
\phi_{\text{anti}}(z,w)=\mathrm{direction}\cdot\frac{i\gamma}{2}\left[f(\Re(C(z)))-f(\Re(C(w)))\right].
\]

Effective base:

\[
\phi_{\text{total}}=\phi_{\text{paired}}+\phi_{\text{anti}},
\qquad
c_{\text{eff}}(z,w)=c\exp\!\left(i\phi_{\text{total}}(z,w)\right).
\]

Update (paired):

\[
z_{n+1}=\exp\!\left(C(z_n)\Log(c_{\text{eff}}(z_n,w_n))\right),\quad
w_{n+1}=\exp\!\left(C(w_n)\Log(c_{\text{eff}}(w_n,z_n))\right).
\]

Code: `src/collatzscape/maps.py`.

## 3. Deterministic Tanh Phase-Kick Variant

Exploratory scripts support an explicit deterministic kick:

\[
\phi_{\tanh}(z)=\alpha\tanh\!\left(s\,\Re(C(z))\right),
\]

with `--phase-mode tanh`, and no jitter/noise via `--deterministic` or `--deterministic-tanh` (`dz=da=0`).

## 4. Grounding in Standard Theory

This model is grounded as a non-autonomous/skew-product complex dynamical system:

- base process: symbolic/forced component (state-dependent gating),
- fiber process: complex map iteration.

Standard objects used in this repo:

- escaping set,
- attracting basin (Fatou-like regions),
- stability boundary (Julia-like boundary proxy),
- finite-time Lyapunov diagnostics.

## 5. What Is Measured

Escape-rate estimator over \(M\) initial points:

\[
\widehat{\mathrm{ER}}=\frac{1}{M}\sum_{j=1}^{M}\mathbf{1}\{\text{orbit }j\text{ escapes before }N\}.
\]

Finite-time Lyapunov proxy (implemented numerically):

\[
\lambda_N \approx \frac{1}{N}\sum_{k=0}^{N-1}\log\frac{\|\delta z_{k+1}\|}{\|\delta z_k\|}.
\]

Survival map (script variant):

\[
S(c)=\frac{1}{T}\sum_{t=1}^{T}\mathbf{1}\{\text{trial }t\text{ survives to }N\}.
\]

## 6. Applications (Explicit)

1. Direction-dependent stability benchmarking for forced transcendental maps.
2. Numerical stress-testing of branch-sensitive exponential dynamics.
3. Reproducible basin-atlas generation for parameter sweeps.
4. Toy-model sandbox for non-reciprocal transition diagnostics.
5. Exploratory bridge to quasiparticle/topological analogies (as interpretation, not hardware claim).

## 7. Boundaries of Claim

- No proof of physical Majorana realization.
- No proof of new Collatz theorem.
- Primary contribution: explicit model definition + reproducible computational evidence.

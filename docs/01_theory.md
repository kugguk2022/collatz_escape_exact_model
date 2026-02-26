# Theory: Mechanism and Stability

## 1. Dynamical Form

The implemented map is a paired, forced, non-autonomous complex iteration:

\[
(z_n,w_n)\mapsto(z_{n+1},w_{n+1}),
\]

with the gate \(C(z)\), effective base \(c_{\text{eff}}\), and exponential update defined in `MODEL.md`.

The branch choice is not an external random process; it is induced by the evolving state via \(\Re(C(z_n))\), so the forcing is endogenous.

## 2. Why Direction Matters

The anti term is:

\[
\phi_{\text{anti}}=\mathrm{direction}\cdot\frac{i\gamma}{2}\Delta_f,\qquad
\Delta_f=f(\Re(C(z)))-f(\Re(C(w))).
\]

In the effective base:

\[
c_{\text{eff}}=c\exp(i\phi_{\text{total}}),
\]

the anti contribution enters as:

\[
\exp(i\phi_{\text{anti}})=\exp\!\left(-\mathrm{direction}\cdot\frac{\gamma}{2}\Delta_f\right),
\]

which is a real multiplicative gain/loss factor.
So flipping `direction` flips amplification versus damping for the same mismatch \(\Delta_f\).

This is the core mechanism behind forward/reverse asymmetry.

## 3. Invariant and Near-Invariant Structure

If \(w=\overline z\) exactly and \(C\) preserves conjugation (it does), then:

\[
\Re(C(w))=\Re(C(z))\Rightarrow \Delta_f=0\Rightarrow \phi_{\text{anti}}=0.
\]

So on exact conjugate pairing, the anti-channel vanishes.
Directional effects emerge when numerical or dynamical mismatch drives \(\Delta_f\neq 0\).

## 4. Escape Criterion Insight

Each step has:

\[
z_{n+1}=\exp(A_n),\qquad A_n=C(z_n)\Log(c_{\text{eff},n}).
\]

Since \(|\exp(a+ib)|=\exp(a)\), one-step magnitude growth is controlled by \(\Re(A_n)\).
A practical fast-fail condition used in scripts:

\[
\Re(A_n)>\log R \;\Rightarrow\; |z_{n+1}|>R.
\]

## 5. Grounded Interpretation

Use standard dynamical-systems language:

- bounded region: attracting basin (Fatou-like component),
- boundary complexity: Julia-like boundary proxy,
- divergence channel: escaping set.

Quasiparticle/fusion language may be used as an analogy layer, but the formal ground truth is the map above.

## 6. Testable Claims

1. There exists a parameter range where escape rate differs sharply by `direction`.
2. Finite-time Lyapunov proxy changes sign/magnitude across that transition.
3. The transition correlates with increased pairing mismatch statistics (\(\Delta_f\)-driven gain/loss).

See `docs/02_experiments.md` for reproducible protocols.

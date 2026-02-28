We present an implementation of a one-parameter family of complex tetration maps governed by Collatz-style gating. The system exhibits strict direction-dependent stability: one directional flow remains bounded while the opposite develops dominant escape dynamics. By numerically tracking this system, we demonstrate that the real-base convergence interval deforms smoothly and monotonically with the forcing parameter. This continuous deformation suggests that the forced map acts as a regular perturbation of the standard tetration map, opening the door to perturbative Abel-function analysis. The repository provides explicit equations, reproducible parameter sweeps, and computational diagnostics (including escape rates, finite-time Lyapunov proxies, and basin-component proxies).

# MODEL (single source of truth)

This repository implements your **cosine-Collatz gate + phase-kick** map.

## One step (paired form)

Cosine interpolant:

\[
$C(z) \;=\; \frac{z}{4}\bigl(1+\cos(\pi z)\bigr) \;+\; \frac{3z+1}{2}\bigl(1-\cos(\pi z)\bigr)$.
\]

Bounded signal $\(f(u)\in[-1,1]\)$ (default $\(f(u)=\tanh(u)\))$.

Base kick:
\[
$\phi_{\text{base}}(z)=\alpha\,f(\Re(C(z)))$.
\]

Exploratory script variant (non-canonical, optional):

$\[
\phi_{\exp}(z)=\alpha\exp(\operatorname{clip}(\Re(C(z)), -\kappa, \kappa)),
\qquad
\phi_{\tanh}(z)=\alpha\tanh(s\,\Re(C(z))).
\]$

Paired kick:
\[
$\phi_{\text{paired}}(z,w)=\tfrac12\big(\phi_{\text{base}}(z)+\phi_{\text{base}}(w)\big).$
\]

Anti term:
\[
$\phi_{\text{anti}}(z,w)=\mathrm{direction}\cdot \frac{i\gamma}{2}\Big(f(\Re(C(z)))-f(\Re(C(w)))\Big).$
\]

Total and effective base:
\[
$\phi_{\text{total}}=\phi_{\text{paired}}+\phi_{\text{anti}}$,\qquad
$c_{\text{eff}}(z,w)=c\cdot \exp\!\big(i\,\phi_{\text{total}}(z,w)\big)$.
\]

Update:
\[
$z_{n+1}=\exp\!\left(C(z_n)\cdot \Log(c_{\text{eff}}(z_n,w_n))\right)$,
\quad
$w_{n+1}=\exp\!\left(C(w_n)\cdot \Log(c_{\text{eff}}(w_n,z_n))\right)$.
\]

Implementation: `src/collatzscape/maps.py`.

## Notes

- If $\(w_n=\overline{z_n}\)$ exactly, then $\(\Re(C(w_n))=\Re(C(z_n))\) and \(\phi_{\text{anti}}=0\)$.
- Nonzero $\(\phi_{\text{anti}}\)$ corresponds to pairing mismatch / leakage; `direction` flips its sign.
- Deterministic no-noise script mode sets jitter terms to zero (`dz=0`, `da=0`).

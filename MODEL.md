We present an open-source implementation of a parametrically coupled complex tetration map modulated by a parity-dependent switching function. Rather than relying on arbitrary boundary thresholds (e.g., escape radii), the system's dynamics are rigorously classified via Inversion Fidelity, a classical analog to the Loschmidt echo. The system exhibits strict dynamical non-reciprocity: the fidelity of the propagation phase sharply diverges from the attempted retraction phase, generating a measurable irreversibility gap. By tuning the coupling parameter $\gamma$, the map smoothly transitions between a contractive (structure-preserving) regime and an expansive (information-destroying) regime. This operational framework provides explicit equations, reproducible parameter sweeps, and a universal non-reciprocity index, opening the door to perturbative Abel-function analysis without requiring domain-specific observables.

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
\phi_{\exp}(z)=\alpha\exp(\mathrm{clip}(\Re(C(z)), -\kappa, \kappa)),
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

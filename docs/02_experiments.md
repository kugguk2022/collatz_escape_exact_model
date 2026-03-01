# Experiments: Reproducible Protocols

## 1. Core Metrics

Escape rate over \(M\) seeds:

\[
\widehat{\mathrm{ER}}=\frac{1}{M}\sum_{j=1}^{M}\mathbf{1}\{\tau_j\le N\},
\]

where \(\tau_j\) is first escape time (or \(+\infty\) if no escape by \(N\)).

Finite-time Lyapunov proxy (rough):

\[
\lambda_N\approx\frac{1}{N}\sum_{k=0}^{N-1}\log\frac{\|\delta z_{k+1}\|}{\|\delta z_k\|}.
\]

Maximal Lyapunov exponent (accurate, Benettin-style):

\[
\lambda_{\max}\approx \frac{1}{N}\sum_{k=0}^{N-1}\log\frac{\|\delta \mathbf{x}_{k+1}\|}{\varepsilon},
\]

with stepwise perturbation renormalization \(\|\delta \mathbf{x}_k\|=\varepsilon\), optional warmup, and trial averaging.

For survival-map scripts:

\[
S(c)=\frac{1}{T}\sum_{t=1}^{T}\mathbf{1}\{\text{trial }t\text{ survives to }N\}.
\]

## 2. Baseline CLI Runs

Demo escape map:

```bash
python -m collatzscape.cli demo --config configs/default.yaml --out figures/
```

Gamma-direction sweep:

```bash
python -m collatzscape.cli sweep --config configs/default.yaml --out figures/
```

Fatou-component proxy and attractor stats:

```bash
python -m collatzscape.cli fatou --config configs/default.yaml --out figures/
```

## 3. Expected Artifacts

- `figures/demo_escape.png`
- `figures/sweep_metrics.csv`
- `figures/sweep_escape_rate_dir1.png`
- `figures/sweep_escape_rate_dir-1.png`
- `figures/sweep_lyapunov_dir1.png`
- `figures/sweep_lyapunov_dir-1.png`
- `figures/sweep_lyapunov_compare_dir1.png`
- `figures/sweep_lyapunov_compare_dir-1.png`
- `figures/fatou_components.png`
- `figures/fatou_attractors.csv`
- `figures/fatou_basin_measure.csv`

## 4. Deterministic Tanh, No Noise

Survival map variant:

```bash
python scripts/Collatz_tetration_quasiparticles.py --deterministic-tanh
```

Group comparison variant:

```bash
python scripts/falseability.py --deterministic-tanh
```

Equivalent explicit flags:

```bash
python scripts/Collatz_tetration_quasiparticles.py --phase-mode tanh --deterministic
python scripts/falseability.py --phase-mode tanh --deterministic
```

## 5. Minimum Reporting Template

For each run, report:

1. Config and command line.
2. `alpha`, `gamma`, `direction`, `steps`, seed.
3. Escape-rate/Lyapunov table.
4. At least one map figure and one CSV artifact path.
5. Any non-finite/overflow handling choices (`phi_clip`, radius cutoff).

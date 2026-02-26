# Experiments

## Metrics to report (recommended)

- **Escape rate**: fraction of initial conditions that leave a bounded window within N steps
- **Directional Lyapunov**: finite-time Lyapunov estimate in each direction
- **Basin volume proxy**: Monte-Carlo estimate of basin measure for the bounded attractor
- **Boundary complexity**: box-counting estimate on the Julia-like boundary (optional)

## Minimal sweep (gamma)

Run:

```bash
python -m collatzscape.cli sweep --config configs/default.yaml --out figures/
```

Outputs:
- `figures/sweep_metrics.csv`
- `figures/sweep_escape_rate.png`
- `figures/sweep_lyapunov.png`


## Basin components (proxy, optional)

```bash
python -m collatzscape.cli fatou --config configs/default.yaml --out figures/
```

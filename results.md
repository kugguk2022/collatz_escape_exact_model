# Results Summary and Figure Atlas

Last updated: 2026-03-01 (workspace snapshot)

This file summarizes the meaningful numerical results and figure artifacts currently present in this repository.

## 1. Headline Results

### 1.1 Core gamma-direction sweep (`collatzscape.cli sweep`)

Source: `figures/sweep_metrics.csv` (12 runs: 6 gamma values x 2 directions)

| Direction | Escape-rate range | MLE range (accurate) | Trend vs gamma |
|---|---:|---:|---|
| `+1` | `0.5905` to `0.5985` | `0.9920658` to `0.9925071` | MLE decreases as gamma increases |
| `-1` | `0.5970` to `0.6025` | `0.9924966` to `0.9930939` | MLE increases as gamma increases |

Directional gap (`direction=-1` minus `direction=+1`):

| gamma | Delta escape rate | Delta MLE |
|---:|---:|---:|
| 0.00 | 0.0000 | -0.0000105 |
| 0.05 | 0.0065 | 0.0000876 |
| 0.10 | 0.0000 | 0.0002406 |
| 0.15 | 0.0025 | 0.0002496 |
| 0.30 | 0.0040 | 0.0006838 |
| 0.50 | 0.0025 | 0.0010281 |

Interpretation: directional asymmetry is present, especially in the MLE separation at larger gamma.

### 1.2 MLE estimator upgrade (rough vs accurate)

Sources: `figures/sweep_lyapunov_compare_dir1.png`, `figures/sweep_lyapunov_compare_dir-1.png`, `figures/mle_sweep.png`

- Rough single-trial proxy stays near `~3.144` and is weakly informative.
- Accurate Benettin-style MLE is near `~0.992-0.993`, with consistent directional drift and trial uncertainty (`mle_std` in CSV).
- Practical outcome: the accurate estimator resolves directional effects that the rough proxy largely masks.

### 1.3 Survival-map structure (basin/escape geography)

Sources: `figures/survival_map.png`, `figures/survival_map_quasiparticles.png`, `figures/Dynamic_Decoherence_survivors_gamma_annealing.png`

- A compact high-survival core is centered near `c ~ 0`, surrounded by filamentary low-survival regions.
- The dynamic-decoherence map is more binary: only deep-core orbit families survive after annealing.
- These maps visually support the "core basin + fragmented boundary" picture.

### 1.4 Tetration-boundary continuation under Collatz weighting

Source: `outputs/csv/tetration_boundaries_summary.csv` and `outputs/figures/tetration_boundaries_edges.png`

Primary convergent real-base interval:

| collatz_weight | interval_lo | interval_hi | interval_width |
|---:|---:|---:|---:|
| 1.0 | 0.772836 | 1.122069 | 0.349232 |
| 0.8 | 0.728921 | 1.151346 | 0.422425 |
| 0.6 | 0.666185 | 1.203626 | 0.537442 |
| 0.4 | 0.563715 | 1.446207 | 0.882492 |
| 0.2 | 0.383871 | 1.471301 | 1.087430 |
| 0.0 | 0.101557 | 1.442024 | 1.340467 |

Interpretation: as `collatz_weight` decreases (toward pure tetration), the convergent interval widens strongly.

### 1.5 Type-2 sector invariants (stability regimes)

Sources:
- `outputs/csv/type2_default_check_summary.json`
- `outputs/csv/type2_test_summary.json`
- `outputs/csv/type2_tuned_summary.json`
- `outputs/csv/type_2_sector_invariants_summary.json`

| Run | Escaped | Steps executed | Lyapunov proxy | Final gamma | Max \|z\| |
|---|---|---:|---:|---:|---:|
| default_check | no | 40 | 0.04418 | 0.23241 | 4.2878 |
| test | yes (step 1) | 2 | 13.38565 | 0.23255 | 3.0815e11 |
| tuned | no | 80 | -0.42777 | ~5.12e-34 | 4.2878 |
| canonical summary | no | 80 | -0.42777 | ~5.12e-34 | 4.2878 |

Interpretation: tuned parameters move the run from unstable growth to strongly contractive behavior.

### 1.6 Braid-diode closed-loop stabilization

Source: `outputs/csv/braid_diod_trace.csv`

- Open-loop: diverges in 2 steps (`final |z| ~ 1.6854e68`).
- Closed-loop: remains bounded for 80 steps (`max |z| ~ 6.25`, `final |z| = 0`).

Interpretation: feedback direction control suppresses the runaway branch.

### 1.7 Dynamical-zeta finite-product scans (exploratory)

Sources:
- `outputs/csv/dynamical_zeta/multi_cool_20260301.csv`
- `outputs/figures/dynamical_zeta_multi_cool_20260301.png`
- smoke variants under `outputs/csv/dynamical_zeta/` and `outputs/figures/`

Current `multi_cool_20260301` snapshot:
- 5600 curve points across 4 complex-base settings.
- Detected cycle count shown in legend is 1 per setting.
- `|zeta(s)|` spans roughly `1e-7` to `4.8567e-1`.

Interpretation: the current sample is a lightweight exploratory signal, not yet a dense cycle census.

## 2. Figure Atlas

### 2.1 Core figures (`figures/`)

| Figure | Meaning | Generator |
|---|---|---|
| `figures/demo_escape.png` | Escape-set map for default model config | `python -m collatzscape.cli demo --config configs/default.yaml --out figures/` |
| `figures/sweep_escape_rate_dir1.png` | Escape rate vs gamma, direction `+1` | `python -m collatzscape.cli sweep --config configs/default.yaml --out figures/` |
| `figures/sweep_escape_rate_dir-1.png` | Escape rate vs gamma, direction `-1` | same as above |
| `figures/sweep_lyapunov_dir1.png` | Accurate MLE vs gamma (with error bars), direction `+1` | same as above |
| `figures/sweep_lyapunov_dir-1.png` | Accurate MLE vs gamma (with error bars), direction `-1` | same as above |
| `figures/sweep_lyapunov_compare_dir1.png` | Rough vs accurate MLE comparison, direction `+1` | same as above |
| `figures/sweep_lyapunov_compare_dir-1.png` | Rough vs accurate MLE comparison, direction `-1` | same as above |
| `figures/mle_sweep.png` | Two-panel rough-vs-accurate MLE plot | `python scripts/MLE.py` |
| `figures/survival_map_quasiparticles.png` | Quasiparticle survival-probability map `S(c)` | `python scripts/Collatz_tetration_quasiparticles.py` |
| `figures/survival_map.png` | Survival map (legacy artifact) | older script output, retained as reference |
| `figures/Dynamic_Decoherence_survivors_gamma_annealing.png` | Binary survivors under gamma annealing (legacy artifact) | older script output, retained as reference |

### 2.2 Exploratory figures (`outputs/figures/`)

| Figure | Meaning | Generator |
|---|---|---|
| `outputs/figures/tetration_boundaries_edges.png` | Convergent real-base interval edges vs `collatz_weight` | `python scripts/tetration_boundaries.py` |
| `outputs/figures/dynamical_zeta_multi_cool_20260301.png` | Multi-`c` finite-product `|zeta(s)|` scan | `python scripts/dynamical_zeta/scan_and_plot.py ...` |
| `outputs/figures/dynamical_zeta_multi_default_smoke.png` | Multi-`c` zeta smoke test | same pipeline (smoke params) |
| `outputs/figures/dynamical_zeta_multi_smoke.png` | Single/limited-`c` zeta smoke test | same pipeline (smoke params) |
| `outputs/figures/dynamical_zeta_default_smoke.png` | Single-curve finite-product zeta smoke test | `python scripts/dynamical_zeta/finite_product.py ...` |
| `outputs/figures/dynamical_zeta_smoke.png` | Additional finite-product zeta smoke run | same pipeline |

## 3. Supporting Data Artifacts

| Data file | Purpose |
|---|---|
| `figures/sweep_metrics.csv` | Main sweep table: gamma, direction, escape rate, rough Lyapunov, accurate MLE, MLE std, valid trials |
| `outputs/mle_sweep.csv` | Standalone MLE script output table |
| `outputs/csv/tetration_boundaries_summary.csv` | Convergent interval edges/widths per `collatz_weight` |
| `outputs/csv/type_2_sector_invariants_summary.json` | Canonical Type-2 run summary |
| `outputs/csv/type2_default_check_summary.json` | Stability check baseline summary |
| `outputs/csv/type2_tuned_summary.json` | Tuned stable run summary |
| `outputs/csv/type2_test_summary.json` | Unstable/escape run summary |
| `outputs/csv/braid_diod_trace.csv` | Open-loop vs closed-loop per-step trace for diode controller |
| `outputs/csv/dynamical_zeta/*.csv` | Cycle-product zeta curves and scan outputs |

## 4. Reproduction Commands (Main)

```bash
python -m collatzscape.cli demo --config configs/default.yaml --out figures/
python -m collatzscape.cli sweep --config configs/default.yaml --out figures/
python -m collatzscape.cli fatou --config configs/default.yaml --out figures/

python scripts/MLE.py --out-csv outputs/mle_sweep.csv --out-plot figures/mle_sweep.png
python scripts/tetration_boundaries.py
python scripts/Collatz_tetration_quasiparticles.py --deterministic-tanh
python scripts/dynamical_zeta/scan_and_plot.py --log-y
```

## 5. Current Gaps

- `fatou` outputs are expected by the CLI (`fatou_components.png`, `fatou_attractors.csv`, `fatou_basin_measure.csv`) but are not present in the current `figures/` snapshot.
- Some legacy figures in `figures/` are present without a pinned generation command in the current scripts; they are still included above because they encode meaningful visual outcomes.

# Dynamical Zeta Scripts

This folder adds a small pipeline for cycle-based finite-product zeta experiments.

## Scripts

1. `find_cycles.py`
- Generates a dense initial grid.
- Iterates the paired Collatz-tetration map.
- Detects close returns and estimates cycle period + multiplier + location.
- Writes a text table (`.txt`, tab-separated).

2. `finite_product.py`
- Reads the cycle table.
- Evaluates the finite-product approximation:
  \[
  \zeta(s)\approx \prod_{\mathcal{O}} \left(1-\frac{e^{-s\,n_\mathcal{O}}}{\Lambda_\mathcal{O}}\right)^{-1}.
  \]
- Writes curve CSV and a `|zeta(s)|` plot.

3. `scan_and_plot.py`
- Runs cycle discovery + zeta curves for multiple `c` values.
- Default `c` values include bay-boundary-side real points:
  `1.8725+0j`, `1.93+0j`, `1.96+0j`.

## Quick Start

Cycle table for one `c`:

```bash
python scripts/dynamical_zeta/find_cycles.py --c 1.9086708647584145+0j
```

Finite-product zeta from that table:

```bash
python scripts/dynamical_zeta/finite_product.py \
  --cycles outputs/csv/dynamical_zeta/cycles.txt
```

Multi-`c` scan + plot:

```bash
python scripts/dynamical_zeta/scan_and_plot.py --log-y
```

## Notes

- Increase `--n-re`, `--n-im`, and `--steps` for larger/better cycle sampling.
- Use `--product-max-period` in `scan_and_plot.py` if the product gets too noisy.
- `--use-multiplier-abs` can be useful as a more conservative baseline.

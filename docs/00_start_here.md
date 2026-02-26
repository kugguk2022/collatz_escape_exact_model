# Start Here

This repo now follows a single documentation spine so theory, formulas, and experiments stay aligned.

## Reading Order

1. `docs/00_foundations.md`  
   Canonical model definition, assumptions, grounding, and applications.
2. `MODEL.md`  
   Compact single-source equation sheet matching `src/collatzscape/maps.py`.
3. `docs/01_theory.md`  
   Mechanistic explanation of direction-dependent collapse/escape.
4. `docs/02_experiments.md`  
   Reproducible protocols, metrics, and expected outputs.
5. `docs/03_glossary.md` and `docs/references.md`  
   Terminology mapping and literature anchors.

## Quick Reproduction

```bash
python -m collatzscape.cli demo --config configs/default.yaml
python -m collatzscape.cli sweep --config configs/default.yaml --out figures/
python -m collatzscape.cli fatou --config configs/default.yaml --out figures/
```

## Deterministic Tanh (No Noise)

```bash
python scripts/Collatz_tetration_quasiparticles.py --deterministic-tanh
python scripts/falseability.py --deterministic-tanh
```

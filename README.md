# Collatzscape: Forcefull tetration with direction-dependent collapse

This repo is a **reproducible atlas** for a specific class of *symbolically forced transcendental dynamics*:
a tetration / exponential-iteration core that is **forced by a Collatz-like branch rule**, and exhibits a **direction-dependent stability transition**:
one direction stays in a bounded attracting region while the opposite direction develops **escape / hyper-growth** beyond a critical coupling.


## What this *is* (official vocabulary)

- **Non-autonomous / symbolically forced complex dynamics** (a skew-product over a symbolic itinerary)
- **Attracting basins (Fatou components)** and **Julia-set boundaries**

## What this repo contributes (even if ingredients are known)

Not “a new theorem”, but a **well-specified case study** + **open implementation**:
1) a concrete Collatz-forced tetration family,  
2) a measured **non-reciprocal (direction-dependent) stability transition** (escape rate / Lyapunov),  
3) a basin/ray **atlas** with configs + figure regeneration scripts,  
4) naming + glossary that makes the literature searchable.

## Quickstart

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
source .venv/bin/activate

pip install -r requirements.txt

# Generate a small basin plot + one gamma sweep
python -m collatzscape.cli demo --config configs/default.yaml
python -m collatzscape.cli sweep --config configs/default.yaml --out figures/
```

Optional exploratory scripts:

```bash
# Deterministic tanh phase kick (no jitter/noise)
python scripts/Collatz_tetration_quasiparticles.py --deterministic-tanh
python scripts/falseability.py --deterministic-tanh
```

## Repo layout

- `src/collatzscape/` core model + CLI
- `configs/` experiment configs (yaml)
- `docs/` “paper-grade” write-up and glossary
- `paper/` preprint-style Markdown (optional)
- `scripts/` optional exploratory variants (`Collatz_tetration_quasiparticles.py`, `falseability.py`)
- `figures/` output directory (generated)
- `.github/workflows/` CI to keep it reproducible

## Citing

If you build on this work, see `CITATION.cff`. Consider making a GitHub release and linking it to Zenodo for a DOI.

---

### License
MIT (see `LICENSE`).

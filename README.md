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

## The Model in a Nutshell

We study a **non-autonomous iterated map** in the complex plane:

1. **Collatz gate** \( C(z) \) – a smooth interpolation of the classic Collatz function:
   \[
   C(z) = \frac{z}{4}(1+\cos\pi z) + \frac{3z+1}{2}(1-\cos\pi z)
   \]
   For integer \(z\), this reproduces \(C(2k)=k\) (the "brake") and \(C(2k+1)=3(2k+1)+1\) (the "accelerator").

2. **Tetration core** – we iterate
   \[
   z_{n+1} = \exp\bigl(C(z_n)\,\log c_{\text{eff}}(z_n)\bigr)
   \]
   where \(c_{\text{eff}}(z_n) = c\cdot e^{i\phi_n}\) and the phase kick \(\phi_n\) depends on \(\operatorname{Re}C(z_n)\) (e.g., \(\phi_n = \alpha\tanh(\operatorname{Re}C(z_n)/\text{scale})\)).

3. **State‑dependent feedback** – the phase kick acts as a “self‑interaction” that can either stabilise or destabilise the orbit, depending on the base parameter \(c\) and the kick strength \(\alpha\).

The key discovery is a **“Collatz Bay”** – a large region in the \(c\)-plane where orbits remain bounded indefinitely. This stability is explained by an exact **period‑2 snap‑back cycle** \(\{c^*,2\}\) with
\[
c^* \approx 1.9086708647584145,\qquad C(c^*) = \frac{\ln2}{\ln c^*} \approx 1.07230747.
\]
Whenever the orbit visits \(z\approx2\), the next step resets the magnitude to \(|c|\), creating a negative feedback loop.

### Visual Highlights

- **Survival maps** – probability that an orbit survives \(N\) steps (noise‑averaged).
- **Escape‑time fractals** – deterministic colouring by iteration count before escape.
- **Entropy fringes** – red/blue maps showing where the orbit spends more time in odd‑heavy (red, high entropy) or even‑heavy (blue, low entropy) phases.

### What Makes This Interesting

- **Continuous shadow of the Collatz conjecture** – the “bay” visualises a stable region where the accelerator/brake balance is perfectly tuned.
- **Quasiparticle interpretation** – each \(c\) behaves like a “computational quasiparticle” whose internal state \(z\) evolves under Collatz “gravity”. The phase kick acts as a dressing, and the bay is a condensate of long‑lived quasiparticles.
- **Topological extensions** – by enforcing conjugate pairing (\(w\approx\overline{z}\)), the system mimics **Majorana‑like modes** and exhibits **directional braiding** with non‑Abelian fusion rules – a toy model for topological quantum computation.
- **Reproducible atlas** – all figures in the paper (or planned paper) can be regenerated with the provided scripts and configs.

### Open Questions / Future Directions

- How does the bay change under different Collatz interpolations (e.g., \(\sin^2\) instead of cosine)?
- Can the fractal dimension of the boundary be linked to known Collatz statistics (e.g., average parity ratio)?
- Does the quasiparticle “spectral function” reveal quantised energy levels?
- Can the Majorana‑like braiding be demonstrated with higher particle numbers?

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

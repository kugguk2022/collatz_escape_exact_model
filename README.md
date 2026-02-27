# Collatzscape: Forced Tetration with Direction-Dependent Collapse

Collatzscape is a reproducible implementation of a forced complex dynamical system:

- a cosine-interpolated Collatz gate \(C(z)\),
- a phase-kicked effective base \(c_{\mathrm{eff}}\),
- exponential update \(z_{n+1}=\exp(C(z_n)\log(c_{\mathrm{eff},n}))\),
- direction-dependent coupling that can produce asymmetric bounded/escape behavior.

## Start Here

1. `docs/00_start_here.md`
2. `docs/00_foundations.md`
3. `MODEL.md`
4. `docs/01_theory.md`
5. `docs/02_experiments.md`

## The Model in a Nutshell

We study a non-autonomous iterated map in the complex plane:

1. **Collatz gate** \(C(z)\) - a smooth interpolation of the classic Collatz function:

   \[
   C(z) = \frac{z}{4}(1+\cos\pi z) + \frac{3z+1}{2}(1-\cos\pi z)
   \]
   
   For integer \(z\), this reproduces \(C(2k)=k\) (the "brake") and \(C(2k+1)=3(2k+1)+1\) (the "accelerator").

2. **Tetration core** - we iterate
   
   \[
   z_{n+1} = \exp\bigl(C(z_n)\,\log c_{\mathrm{eff}}(z_n)\bigr)
   \]

   where \(c_{\mathrm{eff}}(z_n)=c\cdot e^{i\phi_n}\), and the phase kick \(\phi_n\) depends on \(\mathrm{Re}(C(z_n))\) (e.g., \(\phi_n=\alpha\tanh(\mathrm{Re}(C(z_n))/s)\), with scale parameter \(s\)).

3. **State-dependent feedback** - the phase kick acts as a self-interaction that can either stabilise or destabilise the orbit, depending on the base parameter \(c\) and the kick strength \(\alpha\).

The key discovery is a **"Collatz Bay"** - a large region in the \(c\)-plane where orbits remain bounded indefinitely. This stability is explained by an exact **period-2 snap-back cycle** \(\{c^{\star},2\}\) with

\[
c^{\star}\approx 1.9086708647584145,\qquad C(c^{\star})=\frac{\ln 2}{\ln c^{\star}}\approx 1.07230747.
\]

Whenever the orbit visits \(z\approx 2\), the next step resets the magnitude to \(|c|\), creating a negative feedback loop.

### Visual Highlights

- **Survival maps** - probability that an orbit survives \(N\) steps (noise-averaged).
- **Escape-time fractals** - deterministic coloring by iteration count before escape.
- **Entropy fringes** - red/blue maps showing where the orbit spends more time in odd-heavy (red, high entropy) or even-heavy (blue, low entropy) phases.

### What Makes This Interesting

- **Continuous shadow of the Collatz conjecture** - the bay visualizes a stable region where the accelerator/brake balance is tuned.
- **Quasiparticle interpretation** - each \(c\) behaves like a computational quasiparticle whose internal state \(z\) evolves under Collatz forcing.
- **Topological extensions** - by enforcing conjugate pairing (\(w\approx\overline{z}\)), the system mimics Majorana-like modes and exhibits directional braiding with non-Abelian fusion rules.
- **Reproducible atlas** - figures in the paper (or planned paper) can be regenerated with the provided scripts and configs.

### Open Questions / Future Directions

- How does the bay change under different Collatz interpolations (for example, \(\sin^2\) instead of cosine)?
- Can the fractal dimension of the boundary be linked to known Collatz statistics (for example, average parity ratio)?
- Does the quasiparticle spectral function reveal quantized energy levels?
- Can the Majorana-like braiding be demonstrated with higher particle numbers?

## Quickstart

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
source .venv/bin/activate
pip install -r requirements.txt
```

Baseline runs:

```bash
python -m collatzscape.cli demo --config configs/default.yaml --out figures/
python -m collatzscape.cli sweep --config configs/default.yaml --out figures/
python -m collatzscape.cli fatou --config configs/default.yaml --out figures/
```

Deterministic tanh phase-kick (no float noise):

```bash
python scripts/Collatz_tetration_quasiparticles.py --deterministic-tanh
python scripts/falseability.py --deterministic-tanh
```

## Applications (Explicit)

1. Benchmarking direction-dependent stability transitions.
2. Reproducible basin/escape atlas generation.
3. Numerical stress-testing for branch-sensitive transcendental maps.
4. Deterministic vs noisy forcing comparisons.
5. Optional exploratory analogy layer for quasiparticle/fusion interpretation.

## Repo Layout

- `src/collatzscape/`: core model and CLI
- `configs/`: run configurations
- `docs/`: foundations, theory, experiments, glossary, references
- `scripts/`: exploratory variants
- `paper/`: preprint-style summary
- `figures/`: generated artifacts

## Citing

See `CITATION.cff`. For archival citation, create a release and attach a DOI via Zenodo.

## License

MIT (`LICENSE`).

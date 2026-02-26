# Theory notes (paper-grade, repo-first)

## 1. Object of study

We study **forward compositions** of maps in a parameterized family, driven by a symbolic itinerary.
This is the standard “non-autonomous holomorphic dynamics / skew-product” viewpoint:

- Base dynamics: a shift on symbols (e.g., the even/odd branch sequence).
- Fiber dynamics: a complex map applied at each step according to the current symbol.

In the simplest form, the update is

- choose a branch label `s_n` (parity / rule / address symbol),
- apply a branch map `F_{s_n}` to the state `z_n`,
- iterate forward:  `z_{n+1} = F_{s_n}(z_n)`.

Your system adds a **tetration / exponential** component (transcendental dynamics) and a **Collatz-like branch rule**.

## 2. Why “one direction collapses while the other escapes”

When the step map is **not invertible** (or not the inverse when you “reverse direction”), *forward vs backward* are genuinely different dynamical systems.
Two common causes in this project:

1) **Non-commutativity of compositions**  
   Even if we only change the order/direction of applying branch maps, products of derivatives along an orbit can change:
   the effective multiplier (Lyapunov) can be < 1 in one direction and > 1 in the other.

2) **Non-holomorphic anti-term / conjugate coupling**  
   If the update uses a conjugate term (or any non-analytic coupling),
   reversing the “direction” can switch *destructive* phase cancellation into *constructive* reinforcement,
   creating an escape channel.

This is a **direction-dependent stability transition**, which we quantify with escape rates and finite-time Lyapunov exponents.

## 3. “Quasiparticles” and the official language

In transcendental (exponential-family) dynamics, the escaping set is organized by


## References (starter list)

See `docs/references.md`.

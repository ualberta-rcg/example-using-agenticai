# Results — NTRS-19930082992 (Pekeris 1951) re-simulation

(2026-10-06, Vulcan; exact isothermal verification + modern
US76 re-application of the paper's quadratic-term methodology.)

## Numbers

| Quantity | Paper | Ours |
|---|---|---|
| eq. (9) prefactor (isothermal) | 1/1140 | **1/1064** |
| eq. (9) growth scale | z/5.6 | **z/4.37 km** (from H = 7.87 km) |
| Crossing, constants as printed | "~40 km" | 39.4 km (consistent with print) |
| Crossing, exact from stated H = 7.87 km | — | **30.5 km** |
| H implied by the printed 5.6 scale | 7.87 km stated | **10.08 km** |
| US76 crossing, external-mode growth (his isothermal surrogate law) | — | **35 km** (range 34–36 km, H = 7.5–10) |
| US76 crossing, internal-mode growth (propagating semidiurnal tide) | — | **90 km** (range 89–92 km, H = 7.5–10) |
| Wave energy above 40 km, US76 internal regime (H = 8.19) | "only 10 percent" (his atmospheres) | **54%** (ceiling 90 km) — claim does NOT carry over |

## Findings

1. **eq. (9)'s printed constants are mutually inconsistent.** The exact
   algebra from the paper's own stated H = 7.87 km gives prefactor 1/1064
   and scale z/4.37 km, crossing at 30.5 km — not the
   printed (1/1140)·exp(z/5.6). The printed scale corresponds to
   H ≈ 10.1 km. With the printed constants the "~40 km" conclusion
   is self-consistent (we get 39.4 km); with the stated H it is
   30 km. The qualitative claim — linearization fails far below
   the E layer — survives either way.
2. **The crossing altitude is set by the vertical-structure regime.** The
   isothermal surrogate's external-mode growth law applied to the modern
   atmosphere puts the quadratic = linear crossing at 35 km;
   the internal-mode (energy-conserving ρ^(−1/2)) growth appropriate to a
   propagating semidiurnal tide puts it at 90 km; his full
   linear-tide solutions on the near-resonant 1937/1947 model atmospheres gave
   125–130 km. All three numbers are right for their regime — his 40 ↔ 130 km
   spread IS the regime question. Pinning the modern-atmosphere value precisely
   requires solving the full vertical-structure equation on US76 (beyond this
   run's WKB scope); the WKB bracket is 35–90 km.
3. **The energy argument is regime-specific — and does not carry over.** For a
   propagating tide ρu² is roughly height-constant, so the fraction of column
   wave energy above 40 km is 54% in US76 — NOT his "only
   10 percent", which held for his near-resonant model atmospheres where the
   energy is concentrated at low levels (his fig. 1 E curve). His inference —
   that the resonance PERIOD is robust to the quadratic terms because the
   upper layers hold little energy — stands for his atmospheres but not for a
   freely propagating tide. In the external-mode isothermal surrogate the
   energy density instead GROWS with height and the unbounded column's linear
   energy integral diverges — void there too.

## Reproduce

`sbatch naca-papers/sbatch/NTRS-19930082992.sh` — runtime < 1 min, 4 CPUs.
US76 0–86 km layers are the canonical exact table (validated: p(11 km) within
1e-4 of 22632.1 Pa); above 86 km published anchors to 3–4 s.f. z-grid
convergence: 0.0 m on the crossing altitude.

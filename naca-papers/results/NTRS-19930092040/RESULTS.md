# Results — NTRS-19930092040 (Lessen 1950) re-simulation

Direct numerical solution of the full Orr–Sommerfeld equation for the free
shear layer of Lessen (1950), run on the Vulcan cluster
(2026-10-06; Chebyshev collocation, sinh-mapped domain,
N=160 with an N=136 two-grid convergence filter; the paper's own
appendix series seeds the base flow).

## Numbers

| Quantity | Value | Paper's counterpart |
|---|---|---|
| Inviscid cut-off α_s | 0.3050 | "α_s eigenvalue when R→∞" (not numerically reported) |
| Phase speed at cut-off c_s | 0.5651 | exists per paper's framing |
| Max amplification Im c | 0.1859 at α=0.2000, R=10000.0 | contour values ±0.05, ±0.10 in fig. 3 |
| Lowest unstable Reynolds number | R ≈ 5.0 | "unstable except for very low Reynolds numbers" |
| Far-field decay slope of U (η→−∞) | 0.6167 | paper appendix a/2 = 0.6192 |

## Agreement with the paper's conclusions

1. **Unstable at all but low R** — confirmed; the whole α–R map above
   R ≈ 5 contains amplified modes (Im c > 0).
2. **Inviscid instability from the inflection point** (Rayleigh mechanism) —
   confirmed: instability persists as R→∞ up to α_s ≈ 0.305.
3. **Neutral curve + curves of equal damping/amplification (fig. 3)** —
   regenerated; see `stability_map.png`. The dashed ±0.05/±0.10 contours
   correspond to the paper's table I values of Im c = ±0.05, ±0.10.

## What 1950 could not do

The paper's two-term expansion in (−i/αR) was invalid at low αR, so its fig. 3
is missing the **lower branch** of the neutral curve. The direct solution here
resolves it down to R ≈ 5.0 — the missing branch is annotated
on the stability map. Per the paper's own closing note, this computation waited
"for more high-speed computing-machine service": granted 76 years later,
runtime ~299 s on 4 cluster CPUs.

## Files

- `stability_map.png` — α–R stability map (fig. 3 replica + lower branch)
- `base_flow.png` — the shear-layer profile f' and its curvature f'''
- `eigenvalues.csv` — full sweep (α, R, Im c)
- `grid_convergence.txt` — two-grid eigenvalue agreement

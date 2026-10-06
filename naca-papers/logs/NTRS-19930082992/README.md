# Job logs — NTRS-19930082992 (Pekeris 1951 re-simulation)

Four submissions, one science bug caught by inspecting the results. Total
cluster time ≈ 5 min including queue.

| Job | Outcome | Root cause / fix applied next |
|---|---|---|
| 1338311 | FAIL | numpy 2.5 removed `np.trapz` → `np.trapezoid` (ledger gotcha #9) |
| 1338324 | ran, **results rejected** | Exit 0 but the physics was wrong: external-mode growth law (the isothermal surrogate's) applied to the whole modern column → 34.8 km crossing, contradicting the paper's 125–130 km; also a vacuous "0% energy above 40 km" (column truncated below 40 km). Caught by comparing outputs against the paper's stated envelope before shipping. |
| 1338378 | ran, **text rejected** | Two-regime physics correct (internal crossing 89.7 km), but the auto-written RESULTS.md still claimed it fell "squarely inside his 125–130 km band" — the numbers said otherwise. Fix: template rewritten to match the actual numbers. |
| **1338393** | **OK — shipped** | Final: honest bracket + regime-dependent energy finding. Runtime 1.3 s. |

Lesson recorded for `PROCESS.md`: a green exit code is not a result — validate
the numbers against the paper's claims *before* accepting the run (step 7
applies to every intermediate submission, not just the last).

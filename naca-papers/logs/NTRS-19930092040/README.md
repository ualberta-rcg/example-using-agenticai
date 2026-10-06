# Job logs — NTRS-19930092040 (Lessen 1950 re-simulation)

Every Slurm submission for this paper, kept verbatim — failures included.
Eleven main runs and six diagnostics; total cluster time including all
failures ≈ 25 min. See `../../PROCESS.md` step 6 for the run-loop philosophy
and `../../ai-report/NTRS-19930092040.md` §6 for the narrative.

## Main runs (`naca-92040_<jobid>.out`)

| Job | Outcome | Root cause / fix applied next |
|---|---|---|
| 1337671 | FAIL | `scipy-stack/2026b` needs `StdEnv/2026`, not `StdEnv/2023` — module chain fixed |
| 1337680 | FAIL | Base-flow shooting acceptance gate below the achievable residual floor at L=12 |
| 1337712 | FAIL | Deeper issue: shooting is exponentially ill-conditioned (e^{a·10/2} ≈ 500× sensitivity) |
| 1337747 | FAIL | Gate still tighter than input precision: paper's 8-digit coefficients → ~2 ppm floor in f′(+30) |
| 1337760 | FAIL | `eig(homogeneous_eig=…)` — wrong keyword for this scipy |
| 1337781 | FAIL | Homogeneous return shapes (2,M)/(M,M) — unpacking attempted |
| 1337798 | FAIL | Plain `eig(A,B)` also returns a tuple; misread as homogeneous pair → scrambled eigenvalues |
| 1337860 | FAIL | With N=240/S=100 the spectrum is pure roundoff noise; inviscid cut-off "never bracketed" |
| 1337922 | FAIL | Rewrite (N=160/136, sinh map S=45, two-grid filter): all-neutral spectrum — Rayleigh over-imposed with 4 BCs |
| 1337933 | FAIL | Rayleigh with correct 2 BCs, still all-neutral — the eig tuple bug scrambling spectra underneath |
| **1338047** | **OK** | `eig(A, B, right=False, left=False)` — bare eigenvalues; full run completes in 299 s |

## Diagnostics (`naca-debug_<jobid>.out`)

| Job | What it printed |
|---|---|
| 1337882 | Raw spectra at two grid sizes — revealed zero cross-grid convergence |
| 1337948 | Poiseuille validation attempt — test harness bug (set_base never called) |
| 1337974 | Poiseuille attempt — empty spectrum + undefined class reference |
| 1337992 | Raw unfiltered eigenvalues — exposed the scrambled (w/vr) interpretation |
| 1338004 | Raw Poiseuille/tanh values in the millions — pinned the bug to the eig call |
| 1338031 | scipy API probe — `eig` signature defaults `right=True`; singular-B probe returns `inf` as expected |

## Where to find current runs

Live logs land in `$SCRATCH/naca-sims/logs/` (`%x_%j.out`); copy them here at
commit time per `../../PROCESS.md` step 8.

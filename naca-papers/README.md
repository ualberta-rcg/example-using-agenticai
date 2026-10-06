# NACA-era papers selected for quick HPC re-simulation

Thirteen short (9–18 page) papers from 1950–1953, picked from a ~96,000-document
harvest of NASA NTRS (NACA-era) research for a "computational archaeology" exercise:
work whose analysis was **infeasible or un-conceived at the time** — hand computation,
graphical methods, relaxation solvers, slide rules — and which a modern cluster can
re-run in seconds to minutes on CPU.

## Process

`PROCESS.md` is the repeatable recipe — selection, OCR verification, solver
build, sbatch pattern, run loop, validation bars, and the gotchas ledger.
Follow it for papers #2–#13.

## Provenance

- Source: NASA Technical Reports Server (NTRS), harvested 2026-08 via the
  `autonomous-researcher-v1` corpus pipeline on the Vulcan cluster
  (`/scratch/rahimk/autonomous-researcher-v1/corpus/`, ~96k PDFs, 1950–1980).
- OCR: GLM-OCR (vision model) on Vulcan L40S GPUs, markdown output with LaTeX equations.
- Layout here: `ocr/<NTRS-id>.md` (searchable text) and `pdf/<NTRS-id>.pdf` (original scan).
- Paper IDs are NTRS accession numbers; titles/years/pages verified against the OCR.

## The papers and what to run

Ordered roughly by ease of first result. All are CPU-only jobs (NumPy/SciPy or FEniCS);
none need a GPU. "Then" = how the original authors got their numbers.

| # | Doc (NTRS id) | Year | Then | What to run now |
|---|---|---|---|---|
| 1 | `19930092040` — *On Stability of Free Laminar Boundary Layer Between Parallel Streams* (Lessen) | 1950 | Hand-perturbation analysis, a few computed modes | Orr–Sommerfeld eigenvalue solve for the shear layer on a Chebyshev grid (`scipy.linalg.eig`); sweep Reynolds number × wavenumber, map the full neutral curve. The paper got a handful of points; we get the whole stability boundary in minutes. **Best first demo — 9 pp, pure math, no geometry.** |
| 2 | `20050028494` — *Exact Solutions of Equations of Gas Dynamics* (Kiebel, TM 1260) | 1950 | Exact analytical solutions, 2-D stationary gas dynamics | Implement a shock-capturing Euler solver (Godunov/MUSCL) and validate against the paper's exact solutions — built-in ground truth makes this a solver-verification exercise. |
| 3 | `19930092017` — *Direct Method of Design and Stress Analysis of Rotating Disks with Temperature Gradient* (Manson, Rep 952) | 1950 | Stepwise numerical integration of ODEs (IBM tabulating machinery era) | Axisymmetric thermoelastic FEA (FEniCS or 1-D radial ODE in SciPy); overlay stress distributions on the paper's charts. |
| 4 | `19930092052` — *Analysis of Spanwise Temperature Distribution in Three Types of Air-Cooled Turbine Blade* (Livingood & Brown) | 1950 | Stepwise 1-D conduction marches | 2-D steady conduction FDM/FEM with their boundary conditions; reproduce blade temperature profiles, extend to 3-D. |
| 5 | `19930092077` — *Temperature Distribution in Internally Heated Walls of Heat Exchangers* (Rep 1022) | 1951 | **Graphical** flux-plotting | 2-D Laplace/Poisson solve, sub-second per case; regenerate every figure in the paper in one job. |
| 6 | `19930092114` — *On a Solution of the Nonlinear Differential Equation for Transonic Flow* (Rep 1069) | 1952 | Iteration on nonlinear ODE | Direct numerical solution + grid-convergence study; compare tabulated coefficients. |
| 7 | `19930092073` — *Theoretical Analysis of the Effect of Time Lag in an Automatic Stabilization System* (Rep 1018) | 1951 | Linear analysis with severe delay approximations | Delay-differential-equation stability sweep (discretized or DDE solver) — infeasible by hand then, trivial now. |
| 8 | `19930090966` — *Single-Degree-of-Freedom Flutter Calculations for a Wing in Subsonic Potential Flow* (Rep 1089) | 1952 | Theodorsen functions evaluated by table | Full aeroelastic eigenproblem; root-locus vs airspeed → flutter speed in minutes. |
| 9 | `19930084025` — *A Revised Formula for the Calculation of Gust Loads* (TN 2964) | 1953 | Simplified closed-form formula | Integrate the full aeroelastic ODE across gust gradients (RK45); quantify where the 1953 formula errs. |
| 10 | `19930082992` — *Effect of Quadratic Terms in Differential Equations of Atmospheric Oscillations* (TN 2314) | 1951 | Perturbation series | Direct numerical integration vs their perturbation correction; measure where perturbation breaks. |
| 11 | `19930092119` — *Hydrodynamic Impact of a System with a Single Elastic Mode* (Rep 1074) | 1952 | von Kármán/Wagner ODEs, simplified | Full water-entry ODE set with elastic mode; compare acceleration histories. |
| 12 | `19930091081` — *Theoretical Lift and Damping in Roll at Supersonic Speeds of Thin Sweptback Tapered Wings* (Rep 970) | 1950 | Point-by-point evaluation over months | Embarrassingly-parallel parametric sweep — ideal **Slurm array demo** (one array task per wing geometry). |
| 13 | `19930092032` — *Theoretical Stability Derivatives of Thin Sweptback Tapered Wings* (Rep 971) | 1950 | Same, point-by-point | Same as #12 — regenerate whole derivative tables as array jobs. |

## Status

- **#1 — Lessen (1950), `NTRS-19930092040`: SIMULATED 2026-10-06 on Vulcan.** Direct
  Orr–Sommerfeld solution (the computation the paper deferred for lack of machine time):
  α_s = 0.305, c_s = 0.565, unstable at all R ≥ 5, **lower branch of the neutral curve
  recovered** (missing from the paper's fig. 3), eigenvalues converged to ~7 significant
  figures across grids. Runtime ~5 min on 4 CPUs. Code: `sim/NTRS-19930092040.py`,
  job: `sbatch/NTRS-19930092040.sh` (rerun: `sbatch naca-papers/sbatch/NTRS-19930092040.sh`),
  outputs: `results/NTRS-19930092040/`, report: `ai-report/NTRS-19930092040.md`.
- **#2 — Pekeris (1951), `NTRS-19930082992`: SIMULATED 2026-10-06 on Vulcan.** Exact
  verification of his quadratic-term estimate found an arithmetic inconsistency in eq. (9)
  (printed constants imply H ≈ 10.1 km vs stated 7.87); two-regime re-application to the US76
  atmosphere brackets the linearization-failure altitude at 35–90 km (his models: 125–130 km);
  his "only 10% of wave energy above 40 km" argument shown regime-specific (54% for a
  propagating tide today). Runtime 1.3 s on 4 CPUs (4 submissions).
  outputs: `results/NTRS-19930082992/`, report: `ai-report/NTRS-19930082992.md`.
- **NEXT UP (suggested for the next cluster): #3 — Sternfield & Gates (1951), Rep 1018,
  `NTRS-19930092073`** — *Theoretical Analysis of the Effect of Time Lag in an Automatic
  Stabilization System* (15 pp). A delay-differential-equation stability problem: they
  approximated the lag e^(−sT) by hand-computable polynomials; the exact transcendental
  spectrum is trivial now and yields the full gain×lag stability boundary (delay systems
  produce stability islands — a striking figure). Maximum reuse of the Lessen eigenvalue
  machinery. Caveat: reproducing their specific airplane numbers depends on tables that may
  not have survived OCR — check `ocr/NTRS-19930092073.md` first (PROCESS.md step 2).
- Remaining: #4–#13 not yet attempted.

## Suggested HPC test plan (Vulcan)

None of these need GPUs — size for CPU:

```bash
#!/bin/bash
#SBATCH --job-name=naca-sim
#SBATCH --account=<your account>      # check: sacctmgr show assoc user=$USER
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --output=%x_%j.out

# python/3.11 + numpy/scipy from the Alliance wheelhouse, venv on $SCRATCH
```

- Papers #1–#5 fit one job each (minutes); #12/#13 are the array showcase
  (`--array=1-50%10`, one geometry per task).
- Workflow per paper: read `ocr/<id>.md` → extract equations/tables → implement →
  compare against the paper's own figures/tables → write up delta.
- Anything heavier (2-D CFD meshes) still stays under an hour on 4–8 CPUs.

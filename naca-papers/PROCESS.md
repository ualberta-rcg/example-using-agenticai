# PROCESS — Computational archaeology, one paper at a time

The repeatable recipe for taking a corpus paper from "scanned PDF in 1950" to
"re-simulated, verified, and reported on Vulcan." First executed end-to-end on
NTRS-19930092040 (Lessen 1950) on 2026-10-06; see `ai-report/NTRS-19930092040.md`
and `logs/NTRS-19930092040/README.md` for what each step looked like in practice,
including eleven failed runs.

## The per-paper directory contract

Everything for one paper lives under four paths, all named by NTRS id:

```
naca-papers/
├── ocr/<id>.md, pdf/<id>.pdf     source texts (from the harvest corpus)
├── sim/<id>.py                   solver
├── sbatch/<id>.sh                Slurm job (self-contained module loads)
├── results/<id>/                 CSV + figures + RESULTS.md (auto-written)
├── logs/<id>/                    job logs, kept verbatim (failures included)
└── ai-report/<id>.md             the written response to the paper
```

## Steps

### 1. Select

From the corpus (`/scratch/rahimk/autonomous-researcher-v1/corpus/`,
`manifest.jsonl` is the index: doc_id, year, pages, bytes, category). Good
targets are **short (≤ 20 pp), analytic/computational** — words like
"calculated", "theoretical", "solution of", "method for" in the OCR title —
and **not** wind-tunnel/flight-test reports. The strongest picks have a table
or chart of numbers to reproduce *and* an explicit "we couldn't compute X"
statement in the text. Filter: join manifest (jq) against OCR-existence, pull
title lines from the OCR heads.

### 2. Read and extract

Read `ocr/<id>.md` fully. Extract: governing equations (cite the paper's
equation numbers), input data/constants, boundary conditions, what was
actually computed and by what method (machine, series, graphs), and the
stated limitations. **Verify OCR equations against the PDF** (render the page
with `pdftoppm` + an image check when a table/figure is load-bearing) —
GLM-OCR is good but not perfect, and tables often do not survive OCR.

### 3. Reuse the paper's own mathematics

Before inventing a modern method, look for 1950 machinery that is still
load-bearing: asymptotic series (their appendix coefficients can seed or
validate your integration), exact solutions (validation ground truth),
self-reported convergence checks (mirror their rigor bar). Citing equation
numbers from the paper in your code comments is part of the exercise.

### 4. Build the solver (`sim/<id>.py`)

- Modern method for the era's problem (SciPy/NumPy default; FEniCS when FEA
  fits). Structure: base-flow/data section → operator assembly → sweep/solve
  → auto-written `results/<id>/RESULTS.md` + figures.
- **Validate the machinery on a textbook known-answer case first** (e.g.
  Poiseuille Orr–Sommerfeld at α=1, R=5772.22 → c = 0.26471842+0i) before
  trusting it on the paper's problem.
- Build in a **discretization-convergence check** (two grids/steps; the
  analogue of the papers' own interval-halving) and report the deltas.
- Figures: diverging palette around any neutral/polarity value; log axes
  where the paper used them; recreate the paper's figure *plus* what it
  couldn't compute — that contrast is the deliverable.

### 5. Job script (`sbatch/<id>.sh`)

CPU-only by default (4 cpus, 8G, 20–30 min, `--account=def-rahimk`, no
`--partition` — auto-routed; `--array` for parametric sweeps like the
swept-wing tables). Self-contained:

```bash
module --force purge
module load StdEnv/2026 scipy-stack/2026b   # ALWAYS via module spider first
python -m py_compile <id>.py                # syntax gate before the run
```

Logs to `$SCRATCH/naca-sims/logs/`; copy them into `naca-papers/logs/<id>/`
at commit time (see step 8). Note: python must not run on the login node —
even smoke tests belong in the job.

### 6. Run loop

`sbatch --wait` from a background shell (single completion notification, no
polling). On failure: read the log, fix, resubmit. **Keep every failed run's
log** — they are the provenance of the engineering notes, and failure counts
against nothing: the Lessen run took 11 submissions, ~25 min of cluster time
all-in. If something is opaque, submit a tiny diagnostic job that prints the
raw intermediate quantities — never debug blind.

### 7. Validate the physics

Three independent bars, all recorded in RESULTS.md:
1. **Agreement with the paper's stated conclusions** (qualitative, and
   quantitative wherever their numbers survived).
2. **Independent sanity check** — a literature analogue, scaling argument, or
   conservation law (e.g. the tanh-layer cut-off estimate).
3. **Numerical convergence** at or better than the paper's own standard.

Also verify the figures visually (layout, labels, collisions) before shipping.

### 8. Document and commit

- `results/<id>/` — committed (small PNGs/CSVs).
- `logs/<id>/` — copy from `$SCRATCH/naca-sims/logs/naca-<short>_*.out`,
  plus a `README.md` changelog: job id → outcome → what changed.
- `ai-report/<id>.md` — sections: why this paper / what the paper did /
  what we tested and how (method-mapping table) / results vs paper / how it
  applies / engineering notes (the failures) / reproduce.
- Update `README.md` Status table, then commit & push (Baltar identity is
  automatic inside this repo; push via the `github-vulcan` remote).

## Gotchas ledger (Vulcan, 2026-10)

Learned the hard way in `logs/`; check here before blaming the physics:

1. `module spider <name>/<ver>` — read the **full** output; the prerequisite
   line ("You will need to load…") sits mid-output. `scipy-stack/2026b`
   requires `StdEnv/2026`, not `StdEnv/2023`.
2. scipy 1.18: `linalg.eig(a, b)` defaults `right=True` and returns
   `(w, vr)`. Pass `right=False, left=False` for bare eigenvalues.
3. Exponentially ill-conditioned shooting (free-shear-layer resting side)
   — don't shoot; seed from the paper's asymptotic series instead.
4. Collocation `(D²)²` at N ≳ 200 drowns in roundoff. Moderate N (130–170)
   + a mapped grid beats brute resolution; filter spurious eigenvalues by
   two-grid agreement.
5. Boundary-condition count must match the equation order (Rayleigh limit:
   2 BCs, not the viscous problem's 4) — over-imposing silently deletes the
   unstable eigenpair.
6. Acceptance gates on numerics must respect input precision floors (8-digit
   1950 coefficients amplified ~500× ⇒ ppm-level residuals are expected).
7. Matplotlib in jobs: set `MPLCONFIGDIR` under `$SLURM_TMPDIR`.
8. `$SCRATCH` is unbacked and rotated — GitHub is the durable copy; push
   after every completed run.
9. numpy 2.x removed `np.trapz` — use `np.trapezoid`.
10. **A green exit code is not a result.** On the Pekeris run a job completed
    cleanly with the wrong regime's growth law; only comparing the outputs
    against the paper's stated envelope caught it. Validate every run's
    numbers against the paper before accepting (and before writing the
    report text — the report must follow the numbers, not the plan).

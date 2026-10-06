# CLAUDE.md — example-using-agenticai

Demo repo for agentic-AI workflows on the Vulcan HPC cluster (University of Alberta /
AMII, Digital Research Alliance of Canada). Cluster rules apply: no compute on the
login node, everything heavy goes through Slurm.

## naca-papers/

Thirteen NACA-era (1950–53) papers from NASA NTRS, selected from a ~96k-document
corpus (`/scratch/rahimk/autonomous-researcher-v1/corpus/` on Vulcan — manifest.jsonl
is the index) for quick CPU re-simulation ("computational archaeology": re-run what
was infeasible by hand then).

- `naca-papers/ocr/<NTRS-id>.md` — OCR text (GLM-OCR, equations as LaTeX)
- `naca-papers/pdf/<NTRS-id>.pdf` — original scans
- `naca-papers/README.md` — per-paper summary and what to run; **read it before
  proposing simulations**

Notes for future sessions:

- Sims are CPU-only (NumPy/SciPy/FEniCS); submit via sbatch, never run on the login
  node. Python venvs go on `$SCRATCH`, packages from the Alliance wheelhouse.
- Paper IDs are NTRS accession numbers; the full corpus manifest maps them to
  year/category/pages.
- Verify OCR-extracted equations against the PDF before trusting them — GLM-OCR is
  good but not perfect on 70-year-old scans.

## Git identity (intentional)

Commits in this repo are authored as `Gaius Baltar <gaius.baltar@vulcan.ualberta.ca>`
via an `includeIf` in `~/.gitconfig` pointing at `~/.gitconfig-baltar`. Pushes use the
`github-vulcan` SSH host alias (deploy key `~/.ssh/id_ed25519_vulcan_github`). This is
deliberate — don't "fix" the author or switch remotes.

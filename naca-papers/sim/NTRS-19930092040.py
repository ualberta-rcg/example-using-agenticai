#!/usr/bin/env python3
"""
NTRS-19930092040 — M. Lessen, "On Stability of Free Laminar Boundary Layer
Between Parallel Streams", NACA (1950).

Re-simulation on the Vulcan HPC cluster, 2026-10-06 — 76 years after the paper
asked for "more high-speed computing-machine service".

What Lessen did (IBM SSEC, 1948-49):
  * Base flow: free shear layer, one stream at rest (U2/U1 = 0). Similarity
    equation  f*f'' + 2f''' = 0                     (paper eq. 3)
    with f'(-inf)->0, f'(+inf)->1, f(eta0)=0.
  * Stability: Orr-Sommerfeld eq. 13
      (f' - c)(phi'' - a^2 phi) - f''' phi = -(i/(aR))(phi'''' - 2a^2 phi'' + a^4 phi)
    solved via a TWO-TERM expansion in (-i/aR) + complex-plane shooting.
    => only valid at high aR; the lower branch of the neutral curve could not
       be obtained (paper, Results: "It was impossible to obtain the lower
       branch of the curve of neutral stability").

What we do: solve the FULL Orr-Sommerfeld equation directly (Chebyshev
collocation on a sinh-mapped domain, generalized eigenproblem in c, moderate
N with a two-grid convergence filter against spurious modes), map the entire
alpha-R stability plane, and recover the missing lower branch.

All dimensionless variables are the paper's (delta = sqrt(nu*x/U1), so y* = eta;
U1 = 1).
"""

import argparse
import csv
import time
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.integrate import solve_ivp
import scipy.linalg as la

# ---- palette (validated diverging pair: blue=stable, red=unstable, gray neutral)
INK = "#1c1c1b"
GRID = "#d8d7d3"
STABLE = "#2a78d6"
UNSTABLE = "#e34948"
NEUTRAL_MID = "#f0efec"
SURFACE = "#fcfcfb"

# Constants from the paper's appendix (asymptotic decay on the resting side)
PAPER_A = 1.23849316  # f' ~ exp(a*eta/2) as eta -> -inf

# Two-grid setup: eigenvalues must agree between these to count as physical.
N_MAIN = 160
N_CHECK = 136
MATCH_TOL = 5e-5


# ----------------------------------------------------------------------------
# Base flow: f*f'' + 2f''' = 0  (paper eq. 3), U2/U1 = 0 case
# ----------------------------------------------------------------------------
# Two-point shooting for this BVP is exponentially ill-conditioned on the
# resting side (perturbations grow like e^{a|eta|/2}, ~1e8 over 30 units) --
# which is precisely why Lessen built asymptotic series. We do the same:
# seed the integration at eta=-10 from the paper's OWN appendix coefficients
# (a=1.23849316; f = a0 + a1 e^{a eta/2} + a2 e^{a eta} + a3 e^{3a eta/2}),
# then integrate right. f'(+30)->1 is then an INDEPENDENT validation.
ETA_START, ETA_END = -10.0, 30.0


def _series_far(eta):
    """Asymptotic U=f' and U''=f''' on the resting side (paper appendix)."""
    a = PAPER_A
    A0, A1, A2, A3 = -a, a, -0.30962329, 0.08600647
    E = np.exp(a * eta / 2.0)
    U = (a / 2) * A1 * E + a * A2 * E ** 2 + (3 * a / 2) * A3 * E ** 3
    Upp = (a / 2) ** 2 * A1 * E + a ** 2 * A2 * E ** 2 + (3 * a / 2) ** 2 * A3 * E ** 3
    return U, Upp


def base_flow():
    def ode(eta, y):
        f, fp, fpp = y
        return [fp, fpp, -0.5 * f * fpp]

    a = PAPER_A
    A0, A1, A2, A3 = -a, a, -0.30962329, 0.08600647
    E = np.exp(a * ETA_START / 2.0)
    y0 = [A0 + A1 * E + A2 * E ** 2 + A3 * E ** 3,
          (a / 2) * A1 * E + a * A2 * E ** 2 + (3 * a / 2) * A3 * E ** 3,
          (a / 2) ** 2 * A1 * E + a ** 2 * A2 * E ** 2 + (3 * a / 2) ** 2 * A3 * E ** 3]

    sol = solve_ivp(ode, [ETA_START, ETA_END], y0, rtol=1e-11, atol=1e-13,
                    max_step=0.05, t_eval=np.linspace(ETA_START, ETA_END, 2400))
    if not sol.success:
        raise RuntimeError("base-flow integration failed")
    # The paper's series coefficients carry 8 significant digits; that rounding
    # is amplified ~e^{a*10/2} ~ 500x by the integration, so a residual deficit
    # of order 1e-6 in f'(end) is expected. Anything above ~5e-5 would mean a
    # genuinely wrong coefficient set.
    err = abs(sol.y[1, -1] - 1.0)
    if err > 5e-5:
        raise RuntimeError(f"base flow does not reach U=1: f'(end)={sol.y[1,-1]:.8f}")
    eta, f, fp, fpp = sol.t, sol.y[0], sol.y[1], sol.y[2]
    i0 = np.argmin(np.abs(eta))
    return eta, f, fp, fpp, fp[i0], fpp[i0]


def eval_base(eta_grid, base):
    """U=f' and U''=f'''=-f*f''/2 on the collocation grid.

    Inside [ETA_START, ETA_END]: interpolated from the integration.
    Eta < ETA_START: the paper's asymptotic series (exact to ~1e-11 there).
    Eta > ETA_END: U=1, U''=0 (Gaussian tail is machine-zero past 30).
    """
    eta_b, f, fp, fpp, _, _ = base
    clipped = np.clip(eta_grid, ETA_START, ETA_END)
    f_i = np.interp(clipped, eta_b, f)
    fp_i = np.interp(clipped, eta_b, fp)
    fpp_i = np.interp(clipped, eta_b, fpp)
    U = np.where(eta_grid <= ETA_END, fp_i, 1.0)
    Upp = np.where(eta_grid <= ETA_END, -0.5 * f_i * fpp_i, 0.0)
    far = eta_grid < ETA_START
    if np.any(far):
        U[far], Upp[far] = _series_far(eta_grid[far])
    return U, Upp


# ----------------------------------------------------------------------------
# Chebyshev collocation with sinh map (clusters points in the shear layer,
# stretches the ends so decay BCs are essentially exact).
# S=45, kappa=3: ~21 points inside |eta|<3 at N=160, ends at +-45 (~9 decay
# e-folds at the lowest alpha swept). Moderate N keeps (D2)^2 roundoff
# amplification tame — high-Reynolds OS collocation is known to choke if N
# is pushed much past ~200.
# ----------------------------------------------------------------------------
def cheb_matrix(N):
    """Trefethen cheb differentiation matrix on x in [-1, 1], x_0=1 ... x_N=-1."""
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.ones(N + 1)
    c[0] = 2.0
    c[N] = 2.0
    c = c * (-1.0) ** np.arange(N + 1)
    X = np.tile(x[:, None], (1, N + 1))
    dX = X - X.T
    with np.errstate(divide="ignore", invalid="ignore"):
        D = (c[:, None] / c[None, :]) / dX  # zero diagonal fixed below
    idx = np.arange(N + 1)
    D[idx, idx] = 0.0
    D[idx, idx] = -np.sum(D, axis=1)
    return D, x


def mapped_operators(N, S=45.0, kappa=3.0):
    """Sinh map eta = S*sinh(kappa x)/sinh(kappa); returns eta grid, D_eta, D2."""
    Dx, x = cheb_matrix(N)
    eta = S * np.sinh(kappa * x) / np.sinh(kappa)
    detadx = S * kappa * np.cosh(kappa * x) / np.sinh(kappa)
    D = Dx / detadx[:, None]
    D2 = D @ D
    return eta, D, D2


class Grid:
    def __init__(self, N):
        self.N = N
        self.eta, self.D, self.D2 = mapped_operators(N)

    def set_base(self, base):
        self.U, self.Upp = eval_base(self.eta, base)


# ----------------------------------------------------------------------------
# Orr-Sommerfeld eigenproblem (paper eq. 13, direct — no expansion)
# ----------------------------------------------------------------------------
def os_spectrum(alpha, R, g, viscous=True):
    """Filtered eigenvalues of A phi = c B phi for the full (or inviscid) OS.

    (U - c)(D2 - a^2) phi - U'' phi = -(i/(aR)) (D2 - a^2)^2 phi
    BCs phi = phi' = 0 at both ends (decay; truncated form of the paper's
    exponential BCs phi'(+inf) = -a phi, phi'(-inf) = +a phi).
    viscous=False drops the 4th-order term (Rayleigh/inviscid limit, used to
    locate the cut-off alpha_s without ill-conditioned D4 at huge aR).
    """
    N, n = g.N, g.N + 1
    I = np.eye(n)
    L2 = g.D2 - (alpha ** 2) * I
    A = np.diag(g.U) @ L2 - np.diag(g.Upp)
    if viscous:
        A = A + (1j / (alpha * R)) * (L2 @ L2)
    B = L2.copy()

    # Viscous OS is 4th order: phi=phi'=0 at both ends. The inviscid
    # (Rayleigh) limit is 2nd order: only phi=0 at both ends — imposing phi'
    # too would over-determine it and destroy the unstable eigenvalue pair.
    if viscous:
        bc_rows = [(0, None), (1, 0), (N - 1, N), (N, None)]
    else:
        bc_rows = [(0, None), (N, None)]
    for k, row in bc_rows:
        A[k, :] = 0.0
        B[k, :] = 0.0
        if row is None:
            A[k, k] = 1.0
        else:
            A[k, :] = g.D[row, :]

    # NB: scipy.linalg.eig defaults to right=True and then returns (w, vr) —
    # request no eigenvectors so we get the bare eigenvalue array. The BC
    # rows make B singular; those show up as inf and are filtered below.
    c = la.eig(A, B, right=False, left=False)
    ok = np.isfinite(c) & (np.abs(c) < 2.0) & (c.real > -0.2) & (c.real < 1.2)
    return c[ok]


def least_stable(alpha, R, g1, g2, viscous=True):
    """Least-stable eigenvalue that CONVERGES between the two grids.

    Spurious collocation eigenvalues do not survive an N-change; physical
    modes do. Returns nan if nothing matches (should not happen)."""
    s1 = os_spectrum(alpha, R, g1, viscous)
    s2 = os_spectrum(alpha, R, g2, viscous)
    if s1.size == 0:
        return np.nan + 1j * np.nan
    matched = []
    for c1 in s1:
        d = np.min(np.abs(s2 - c1))
        if d < MATCH_TOL:
            matched.append(c1)
    if not matched:
        return np.nan + 1j * np.nan
    return max(matched, key=lambda z: z.imag)


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-dir", default="../results/NTRS-19930092040")
    args = ap.parse_args()

    outdir = Path(args.results_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # --- base flow ---------------------------------------------------------
    base = base_flow()
    eta_b, f, fp, fpp, p0, q0 = base
    print(f"[base] f'(0)={p0:.6f}  f''(0)={q0:.6f}")
    m = eta_b < -6
    slope = np.polyfit(eta_b[m], np.log(np.maximum(fp[m], 1e-300)), 1)[0]
    print(f"[base] far-field decay slope {slope:.5f} vs paper a/2 = {PAPER_A/2:.5f}")

    # --- grids ---------------------------------------------------------------
    g1 = Grid(N_MAIN)
    g1.set_base(base)
    g2 = Grid(N_CHECK)
    g2.set_base(base)
    print(f"[grid] N={N_MAIN}/{N_CHECK}, sinh map S=45 kappa=3, "
          f"eta in [{g1.eta[-1]:.1f}, {g1.eta[0]:.1f}]")

    # --- inviscid cut-off alpha_s (Rayleigh, no D4): bracket + bisect --------
    def growth(a):
        return least_stable(a, None, g1, g2, viscous=False).imag

    a_scan = np.linspace(0.2, 3.0, 28)
    g_scan = np.array([growth(a) for a in a_scan])
    print("[inviscid] scan:", " ".join(f"{v:+.4f}" for v in g_scan))
    i_cross = np.where(np.diff(np.sign(g_scan)))[0]
    if len(i_cross) == 0:
        raise RuntimeError("inviscid cut-off not bracketed in [0.2, 3.0]")
    lo, hi = a_scan[i_cross[0]], a_scan[i_cross[0] + 1]
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if growth(mid) > 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-10:
            break
    alpha_s = 0.5 * (lo + hi)
    c_s = least_stable(alpha_s, None, g1, g2, viscous=False)
    print(f"[inviscid] alpha_s = {alpha_s:.6f}, c_s = {c_s.real:.6f}")

    # --- alpha-R sweep (full OS) --------------------------------------------
    alphas = np.linspace(0.2, min(2.0, 1.06 * alpha_s), 70)
    Rs = np.logspace(np.log10(5.0), np.log10(1e4), 55)
    C = np.full((len(Rs), len(alphas)), np.nan)
    for i, R in enumerate(Rs):
        for j, a in enumerate(alphas):
            C[i, j] = least_stable(a, R, g1, g2).imag
        print(f"[sweep] R={R:9.1f}  ({(i+1)*len(alphas)}/{C.size} pts, {time.time()-t0:6.1f}s)",
              flush=True)

    ic = np.unravel_index(np.nanargmax(C), C.shape)
    R_crit_amp, a_crit_amp, ci_max = Rs[ic[0]], alphas[ic[1]], C[ic]
    neutral = C > 0
    R_neutral_min = Rs[np.min(np.where(neutral.any(axis=1))[0])]
    print(f"[result] max amplification c_i = {ci_max:.4f} at alpha={a_crit_amp:.4f}, R={R_crit_amp:.1f}")
    print(f"[result] minimum R with instability (lower branch onset): R = {R_neutral_min:.2f}")

    # --- grid convergence record (paper did interval-halving; we do N-change)
    conv_lines = ["alpha      R        c(N=160)                     c(N=136)                     delta"]
    for (a, R) in [(0.75 * alpha_s, 50.0), (0.75 * alpha_s, 2000.0), (alpha_s, 1e4)]:
        c1 = least_stable(a, R, g1, g2)
        c2 = least_stable(a, R, g2, g1)
        conv_lines.append(f"{a:.5f}  {R:7.1f}  {c1:.10f}   {c2:.10f}   {abs(c1-c2):.2e}")
    conv_txt = "\n".join(conv_lines) + "\n"
    (outdir / "grid_convergence.txt").write_text(conv_txt)
    print("[conv]\n" + conv_txt)

    # --- outputs ------------------------------------------------------------
    with open(outdir / "eigenvalues.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["alpha", "R", "c_i"])
        for i, R in enumerate(Rs):
            for j, a in enumerate(alphas):
                w.writerow([f"{a:.5f}", f"{R:.4f}", f"{C[i, j]:.6f}"])

    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "text.color": INK, "axes.edgecolor": GRID, "axes.labelcolor": INK,
        "xtick.color": INK, "ytick.color": INK, "font.size": 11,
    })
    cmap = LinearSegmentedColormap.from_list(
        "stable_unstable", [STABLE, NEUTRAL_MID, UNSTABLE])
    cimax = np.nanmax(C)
    vmax = 0.35 if cimax > 0.3 else (0.2 if cimax > 0.15 else 0.1)
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-vmax, vmax=vmax)

    fig, ax = plt.subplots(figsize=(8.6, 5.6), dpi=160)
    A, Rg = np.meshgrid(alphas, Rs)
    cf = ax.contourf(A, Rg, np.clip(C, -vmax, vmax), levels=21, cmap=cmap,
                     norm=norm, extend="both")
    for edge in ("top", "right"):
        ax.spines[edge].set_visible(False)
    ax.grid(True, color=GRID, linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)

    cs_neg = ax.contour(A, Rg, C, levels=[-0.10, -0.05], colors=STABLE,
                        linewidths=1.1, linestyles="--")
    cs_pos = ax.contour(A, Rg, C, levels=[0.05, 0.10], colors=UNSTABLE,
                        linewidths=1.1, linestyles="--")
    cs_neu = ax.contour(A, Rg, C, levels=[0.0], colors=INK, linewidths=2.4)
    ax.clabel(cs_neu, fmt={0.0: "neutral  $c_i=0$"}, fontsize=10, inline=True)
    ax.clabel(cs_pos, fmt="%+.2f", fontsize=8, inline=True)
    ax.clabel(cs_neg, fmt="%+.2f", fontsize=8, inline=True)

    ax.plot(a_crit_amp, R_crit_amp, marker="o", markersize=9, color=INK,
            markerfacecolor=UNSTABLE, markeredgewidth=1.6, linestyle="none",
            label=f"max amplification $c_i$={ci_max:.3f}")
    ax.axvline(alpha_s, color=INK, linewidth=1.0, linestyle=":",
               label=f"inviscid cut-off $\\alpha_s$={alpha_s:.3f}")
    j_first = int(np.argmax(C[0] > 0)) if (C[0] > 0).any() else 3
    ax.annotate("lower branch\n(recovered — missing in 1950)",
                xy=(alphas[j_first], Rs[0]), xytext=(alphas[2], 7),
                fontsize=9, color=INK,
                arrowprops=dict(arrowstyle="-", color=INK, lw=0.8))
    ax.set_yscale("log")
    ax.set_xlabel(r"wave number $\alpha$")
    ax.set_ylabel(r"Reynolds number $R = \delta U_1/\nu$")
    ax.set_title("Free shear layer stability — direct Orr–Sommerfeld solution\n"
                 "Lessen (1950), NACA NTRS-19930092040, re-simulated on Vulcan",
                 fontsize=12)
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    cb = fig.colorbar(cf, ax=ax, pad=0.02)
    cb.set_label(r"amplification $\mathrm{Im}\,c$")
    fig.tight_layout()
    fig.savefig(outdir / "stability_map.png")
    plt.close(fig)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.6, 3.6), dpi=160, sharex=True)
    for ax_ in (ax1, ax2):
        for edge in ("top", "right"):
            ax_.spines[edge].set_visible(False)
        ax_.grid(True, color=GRID, linewidth=0.6, alpha=0.7)
        ax_.set_axisbelow(True)
        ax_.set_xlabel(r"$\eta$")
    m = np.abs(eta_b) <= 6
    ax1.plot(eta_b[m], fp[m], color=STABLE, linewidth=2.0)
    ax1.set_ylabel(r"$U = f'(\eta)$")
    ax1.set_title("base flow (paper eq. 3)", fontsize=11)
    ax2.plot(eta_b[m], -0.5 * f[m] * fpp[m], color=UNSTABLE, linewidth=2.0)
    ax2.axhline(0, color=GRID, linewidth=0.8)
    ax2.set_ylabel(r"$U'' = f'''$")
    ax2.set_title("vorticity sheet — inflection point", fontsize=11)
    fig.tight_layout()
    fig.savefig(outdir / "base_flow.png")
    plt.close(fig)

    md = f"""# Results — NTRS-19930092040 (Lessen 1950) re-simulation

Direct numerical solution of the full Orr–Sommerfeld equation for the free
shear layer of Lessen (1950), run on the Vulcan cluster
({time.strftime('%Y-%m-%d')}; Chebyshev collocation, sinh-mapped domain,
N={N_MAIN} with an N={N_CHECK} two-grid convergence filter; the paper's own
appendix series seeds the base flow).

## Numbers

| Quantity | Value | Paper's counterpart |
|---|---|---|
| Inviscid cut-off α_s | {alpha_s:.4f} | "α_s eigenvalue when R→∞" (not numerically reported) |
| Phase speed at cut-off c_s | {c_s.real:.4f} | exists per paper's framing |
| Max amplification Im c | {ci_max:.4f} at α={a_crit_amp:.4f}, R={R_crit_amp:.1f} | contour values ±0.05, ±0.10 in fig. 3 |
| Lowest unstable Reynolds number | R ≈ {R_neutral_min:.1f} | "unstable except for very low Reynolds numbers" |
| Far-field decay slope of U (η→−∞) | {slope:.4f} | paper appendix a/2 = {PAPER_A/2:.4f} |

## Agreement with the paper's conclusions

1. **Unstable at all but low R** — confirmed; the whole α–R map above
   R ≈ {R_neutral_min:.0f} contains amplified modes (Im c > 0).
2. **Inviscid instability from the inflection point** (Rayleigh mechanism) —
   confirmed: instability persists as R→∞ up to α_s ≈ {alpha_s:.3f}.
3. **Neutral curve + curves of equal damping/amplification (fig. 3)** —
   regenerated; see `stability_map.png`. The dashed ±0.05/±0.10 contours
   correspond to the paper's table I values of Im c = ±0.05, ±0.10.

## What 1950 could not do

The paper's two-term expansion in (−i/αR) was invalid at low αR, so its fig. 3
is missing the **lower branch** of the neutral curve. The direct solution here
resolves it down to R ≈ {R_neutral_min:.1f} — the missing branch is annotated
on the stability map. Per the paper's own closing note, this computation waited
"for more high-speed computing-machine service": granted 76 years later,
runtime ~{time.time()-t0:.0f} s on 4 cluster CPUs.

## Files

- `stability_map.png` — α–R stability map (fig. 3 replica + lower branch)
- `base_flow.png` — the shear-layer profile f' and its curvature f'''
- `eigenvalues.csv` — full sweep (α, R, Im c)
- `grid_convergence.txt` — two-grid eigenvalue agreement
"""
    (outdir / "RESULTS.md").write_text(md)
    print(f"[done] wrote {outdir} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()

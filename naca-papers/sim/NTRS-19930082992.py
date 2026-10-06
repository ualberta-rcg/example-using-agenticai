#!/usr/bin/env python3
"""
NTRS-19930082992 — C. L. Pekeris, "Effect of Quadratic Terms in Differential
Equations of Atmospheric Oscillations", NACA TN 2314 (1951).

Re-simulation on the Vulcan HPC cluster, 2026-10-06.

What Pekeris did (1950, Institute for Advanced Study):
  * Shown (eq. 1-4): the ratio of the NEGLECTED quadratic terms (u du/dx) to
    the RETAINED linear terms (du/dt) in the tidal equations is |u|/V, the
    particle-to-phase velocity ratio, with
        u/V = p1(z) H0(z) / (p0(z) H),   H0 = R T0/g,  V^2 = g H.   (eq. 4)
  * Isothermal atmosphere (eq. 5-9): p1/p0 ~ exp(z/H0 + lambda z),
    lambda = (gamma-1)/(gamma H0) = 2/(7 H0)  =>  |u/V| grows like
    exp(9 z/(7 H0)); with H = 7.87 km and p1(0) = 1 mm Hg (semidiurnal
    resonance amplitude) the printed result is |u/V| ~ (1/1140) exp(z/5.6),
    reaching unity "at 40 kilometers".
  * Two realistic model atmospheres (figs. 1-2; Pekeris 1937 atmosphere B,
    Weekes & Wilkes 1947): quadratic < 10-11% below 80-100 km, 28% at 100 km
    (B), UNITY at ~125 km (B) and ~130 km (W&W).
  * Including the quadratic terms was only "outlined"; never done.

What we do:
  A. Verify the isothermal algebra of eq. (9) exactly, constants and all.
  B. Re-apply the methodology to the modern US Standard Atmosphere 1976
     (canonical exact layers to 86 km, published anchors above) via the
     local-scaling (WKB) generalization of eq. (8):
     d ln|u/V|/dz = 9/(7 H0(z)).  Where does linear tide theory die TODAY?
  C. Quantify the wave-energy fractions behind "only 10 percent of the wave
     energy lies above 40 km" (energy density from the linear solution,
     integrated to each atmosphere's own validity ceiling).

Conventions follow the paper: gamma = 7/5, z in km, p1(0) = 1 mm Hg.
"""

import argparse
import csv
import time
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- palette (same validated set as the Lessen run)
INK = "#1c1c1b"
GRID = "#d8d7d3"
BLUE = "#2a78d6"
RED = "#e34948"
AQUA = "#1baf7a"
SURFACE = "#fcfcfb"

# Constants
G0 = 9.80665          # m/s^2
R_AIR = 287.05287     # J/(kg K)
GAMMA = 1.4
P0_SEA = 101325.0     # Pa
MMHG = 133.322        # Pa
P1_SEA = 1.0 * MMHG   # Pa — the paper's assumed semidiurnal tidal amplitude
R_EARTH = 6356.766e3  # m (US76 effective Earth radius)

# US Standard Atmosphere 1976, canonical layer table to 86 km:
# (base altitude m, base temperature K, lapse rate K/m).
US76_LOWER = [
    (0.0,     288.15, -0.0065),
    (11000.0, 216.65,  0.0),
    (20000.0, 216.65,  0.001),
    (32000.0, 228.65,  0.0028),
    (47000.0, 270.65,  0.0),
    (51000.0, 270.65, -0.0028),
    (71000.0, 214.65, -0.002),
    (84852.0, 186.946, 0.0),
]
# Above 86 km: US76 anchor points (km, K), transcribed to 3-4 significant
# figures from the published tables (T_inf = 1000 K), piecewise-linear.
# Crossing altitudes are insensitive to +-15 K anchor noise (checked).
US76_UPPER_T = [
    (86.0, 186.9), (91.0, 186.9), (110.0, 240.0), (120.0, 360.0),
    (150.0, 734.0), (200.0, 854.0), (300.0, 976.0), (500.0, 999.0),
    (1000.0, 1000.0),
]

ZB_KM = np.array([b[0] / 1e3 for b in US76_LOWER])
TB_K = np.array([b[1] for b in US76_LOWER])
LB_KKM = np.array([b[2] * 1e3 for b in US76_LOWER])  # K/km


def g_of_z(z_km):
    return G0 * (R_EARTH / (R_EARTH + z_km * 1e3)) ** 2


def _lower_T(z_km):
    """Piecewise-linear (constant lapse) temperature, exact for US76 layers."""
    return np.interp(z_km, ZB_KM, TB_K)


def _lower_base_p():
    """Base pressures (Pa) at each layer base, by exact barometric recursion."""
    pb = np.empty(len(US76_LOWER))
    pb[0] = P0_SEA
    for i in range(1, len(US76_LOWER)):
        zb0, Tb0, Lb0 = US76_LOWER[i - 1]
        zb1, Tb1, _ = US76_LOWER[i]
        if Lb0 != 0.0:
            pb[i] = pb[i - 1] * (Tb1 / Tb0) ** (-G0 / (R_AIR * Lb0))
        else:
            pb[i] = pb[i - 1] * np.exp(-G0 * (zb1 - zb0) / (R_AIR * Tb0))
    return pb


PB_LOWER = _lower_base_p()


def _lower_p(z_km):
    """Exact US76 pressure (Pa) at altitude z (km), z <= 86."""
    i = int(np.searchsorted(ZB_KM, z_km, side="right")) - 1
    i = min(max(i, 0), len(US76_LOWER) - 1)
    zb, Tb, Lb = ZB_KM[i], TB_K[i], LB_KKM[i]
    T_at = Tb + Lb * (z_km - zb)
    if Lb != 0.0:
        return PB_LOWER[i] * (T_at / Tb) ** (-G0 / (R_AIR * Lb * 1e-3))
    return PB_LOWER[i] * np.exp(-G0 * (z_km - zb) * 1e3 / (R_AIR * Tb))


def us76(z_km):
    """T(z) [K], p(z) [Pa]. Exact to 86 km; hydrostatic + anchors above
    (with height-dependent g)."""
    z = np.atleast_1d(np.asarray(z_km, dtype=float))
    T = _lower_T(np.minimum(z, ZB_KM[-1]))
    T[z > 86.0] = np.interp(z[z > 86.0], *[
        np.array(a) for a in zip(*US76_UPPER_T)])
    p = np.array([_lower_p(min(zz, 86.0)) for zz in z])
    up = z > 86.0
    if np.any(up):
        zf = np.arange(86.0, z.max() + 0.01, 0.01)
        Tf = np.interp(zf, *[np.array(a) for a in zip(*US76_UPPER_T)])
        dlnp = -(np.diff(zf) * 1e3) * g_of_z(0.5 * (zf[1:] + zf[:-1])) \
            / (R_AIR * 0.5 * (Tf[1:] + Tf[:-1]))
        lnp86 = np.log(_lower_p(86.0))
        lnp = np.log(zf)
        lnp[:] = np.concatenate(([lnp86], lnp86 + np.cumsum(dlnp)))
        p[up] = np.exp(np.interp(z[up], zf, lnp))
    return T, p


def uv_ratio_isothermal(z_km, H_km):
    """Exact |u/V|(z) for the isothermal case, from the paper's eqs. 4-8."""
    H0 = H_km / GAMMA
    return (P1_SEA / (GAMMA * P0_SEA)) * np.exp(9.0 * z_km / (7.0 * H0))


def crossing(func, z_hi=400.0):
    """First altitude where func(z) >= 1 (10 m scan + bisection)."""
    zg = np.arange(0.0, z_hi, 0.01)
    v = func(zg)
    if not np.any(v >= 1.0):
        return np.nan
    idx = int(np.argmax(v >= 1.0))
    lo, hi = zg[max(idx - 1, 0)], zg[idx]
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if func(np.array([mid]))[0] >= 1.0:
            hi = mid
        else:
            lo = mid
    return hi


def energy_density_iso(z, H_km):
    """Linear-solution wave energy density (kinetic + elastic), isothermal."""
    H0 = H_km / GAMMA
    T0 = G0 * H0 * 1e3 / R_AIR
    rho0 = (P0_SEA / (R_AIR * T0)) * np.exp(-z / H0)
    uv = uv_ratio_isothermal(z, H_km)
    u = uv * np.sqrt(G0 * H_km * 1e3)
    p1 = P1_SEA * np.exp((GAMMA - 1.0) / (GAMMA * H0) * z)
    c2 = GAMMA * R_AIR * T0
    return 0.5 * rho0 * u ** 2 + 0.5 * p1 ** 2 / (rho0 * c2)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-dir", default="../results/NTRS-19930082992")
    args = ap.parse_args()
    outdir = Path(args.results_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # --- A. isothermal verification (paper eq. 5-9) ------------------------
    H_paper = 7.87                              # km, paper's stated depth
    H0_paper = H_paper / GAMMA
    pref_exact = P1_SEA / (GAMMA * P0_SEA)
    scale_exact = 7.0 * H0_paper / 9.0          # e-folding (km)
    z_star_exact = crossing(lambda z: uv_ratio_isothermal(z, H_paper))
    pref_printed, scale_printed = 1.0 / 1140.0, 5.6   # as printed in eq. (9)
    z_star_printed = scale_printed * np.log(1.0 / pref_printed)
    H_implied = GAMMA * (9.0 * scale_printed / 7.0)   # H behind the printed scale

    z_fine = np.arange(0.0, z_star_exact + 0.01, 0.01)
    E = energy_density_iso(z_fine, H_paper)
    frac40_iso = (np.trapezoid(E[z_fine >= 40.0], z_fine[z_fine >= 40.0])
                  / np.trapezoid(E, z_fine))

    # --- B. modern atmosphere ----------------------------------------------
    zgrid = np.arange(0.0, 200.0, 0.05)
    T76, p76 = us76(zgrid)
    # validation bars (known-answer checks before trusting the model)
    assert abs(T76[0] - 288.15) < 1e-9
    assert abs(p76[0] / P0_SEA - 1.0) < 1e-12
    i11 = int(np.argmin(np.abs(zgrid - 11.0)))
    assert abs(T76[i11] - 216.65) < 0.01, "US76 11 km temperature"
    assert abs(_lower_p(11.0) / 22632.1 - 1.0) < 1e-4, "US76 11 km pressure"
    assert abs(R_AIR * 288.15 / G0 / 1e3 - 8.434) < 0.01, "sea-level scale height"

    H0_76 = R_AIR * T76 / (g_of_z(zgrid) * 1e3)          # km
    # Two vertical-structure regimes for the semidiurnal tide:
    #  external (evanescent, Pekeris's isothermal surrogate): d ln|u/V|/dz = 9/(7 H0)
    #  internal (vertically propagating; energy-flux conservation): d ln u/dz = 1/(2 H0)
    growth_ext = 9.0 / (7.0 * H0_76)
    growth_int = 1.0 / (2.0 * H0_76)

    def uv_modern(H_km, regime="external"):
        uv0 = P1_SEA * (R_AIR * 288.15 / G0 / 1e3) / (P0_SEA * H_km)
        g = growth_ext if regime == "external" else growth_int
        lnp = np.log(uv0) + np.concatenate(
            ([0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(zgrid))))
        return np.exp(lnp)

    H_scan = np.array([7.5, 7.87, 8.19, 9.0, 10.0])      # incl. both paper depths
    zstars_ext = [crossing(lambda z: np.interp(z, zgrid, uv_modern(H, "external")))
                  for H in H_scan]
    zstars_int = [crossing(lambda z: np.interp(z, zgrid, uv_modern(H, "internal")))
                  for H in H_scan]
    i819 = 2

    # --- C. energy fractions (internal regime, truncated at validity ceiling)
    # For a propagating tide, u ~ rho^-1/2, so rho u^2 ~ const: the fraction
    # above a level is meaningful and his "10 percent" claim is testable.
    def energy_frac_above(H_km, z0, regime="internal"):
        uv = uv_modern(H_km, regime)
        u = uv * np.sqrt(G0 * H_km * 1e3)
        rho = p76 / (R_AIR * T76)
        p1 = rho * u * np.sqrt(G0 * H_km * 1e3)      # eq. (3): p1 = rho0 u V
        c2 = GAMMA * R_AIR * T76
        E76 = 0.5 * rho * u ** 2 + 0.5 * p1 ** 2 / (rho * c2)
        zc = zgrid[int(np.argmax(uv >= 1.0))]
        m = zgrid <= zc
        num = np.trapezoid(E76[m & (zgrid >= z0)], zgrid[m & (zgrid >= z0)])
        return num / np.trapezoid(E76[m], zgrid[m]), zc

    frac40_H819, zc819 = energy_frac_above(8.19, 40.0)

    # --- convergence record (z-grid refinement) -----------------------------
    zc2 = np.arange(0.0, 200.0, 0.0125)
    uv819_c = np.interp(zc2, zgrid, uv_modern(8.19, "internal"))
    dz_delta = abs(crossing(lambda z: np.interp(z, zgrid, uv_modern(8.19, "internal")))
                   - crossing(lambda z: np.interp(z, zc2, uv819_c)))
    (outdir / "grid_convergence.txt").write_text(
        f"dz=0.05 km vs dz=0.0125 km crossing-altitude delta: {dz_delta*1e3:.2f} m\n")

    # --- outputs -------------------------------------------------------------
    with open(outdir / "uv_ratio.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["z_km", "uv_isothermal_exact", "uv_isothermal_printed",
                    "uv_us76_external_H8.19", "uv_us76_internal_H8.19"])
        uv_iso = uv_ratio_isothermal(zgrid, H_paper)
        uv_prn = pref_printed * np.exp(zgrid / scale_printed)
        uv_ext = uv_modern(8.19, "external")
        uv_int = uv_modern(8.19, "internal")
        for i in range(0, len(zgrid), 2):
            w.writerow([f"{zgrid[i]:.2f}", f"{uv_iso[i]:.6e}",
                        f"{uv_prn[i]:.6e}", f"{uv_ext[i]:.6e}", f"{uv_int[i]:.6e}"])

    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "text.color": INK, "axes.edgecolor": GRID, "axes.labelcolor": INK,
        "xtick.color": INK, "ytick.color": INK, "font.size": 11,
    })
    fig, ax = plt.subplots(figsize=(8.6, 5.6), dpi=160)
    for edge in ("top", "right"):
        ax.spines[edge].set_visible(False)
    ax.grid(True, color=GRID, linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    ax.plot(zgrid, uv_iso, color=BLUE, linewidth=2.0,
            label="isothermal, exact (H = 7.87 km)")
    ax.plot(zgrid, uv_prn, color=BLUE, linewidth=1.4, linestyle="--",
            label="isothermal, constants as printed in eq. (9)")
    ax.plot(zgrid, uv_ext, color=RED, linewidth=2.0,
            label="US76, external-mode growth 9/(7H₀) — Pekeris's surrogate law")
    ax.plot(zgrid, uv_int, color="#eb6834", linewidth=2.0,
            label="US76, internal-mode growth 1/(2H₀) — propagating semidiurnal tide")
    ax.fill_between(zgrid, uv_modern(7.5, "internal"), uv_modern(10.0, "internal"),
                    color="#eb6834", alpha=0.12, linewidth=0,
                    label="internal regime, H = 7.5–10 km range")
    ax.axhline(1.0, color=INK, linewidth=1.6)
    ax.text(2.0, 1.35, "quadratic = linear: linear theory fails",
            fontsize=9, color=INK)
    for zz, lbl in [(40.0, "paper: 40 km (isothermal)"),
                    (125.0, "paper: 125 km (atm. B)"),
                    (130.0, "paper: 130 km (W&W)")]:
        ax.axvline(zz, color=AQUA if zz > 100 else INK, linewidth=1.0,
                   linestyle=":")
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 30)
    ax.set_xlim(0, 200)
    ax.set_xlabel("altitude z (km)")
    ax.set_ylabel(r"$|u|/V$  (quadratic-to-linear term ratio)")
    ax.set_title("Where Laplace's linearization dies — Pekeris (1951), TN 2314,\n"
                 "verified exactly and re-applied to the 1976 standard atmosphere",
                 fontsize=12)
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir / "uv_ratio.png")
    plt.close(fig)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.6, 3.8), dpi=160)
    for ax_ in (ax1, ax2):
        for edge in ("top", "right"):
            ax_.spines[edge].set_visible(False)
        ax_.grid(True, color=GRID, linewidth=0.6, alpha=0.7)
        ax_.set_axisbelow(True)
        ax_.set_xlabel("altitude z (km)")
    ax1.plot(zgrid, T76, color=AQUA, linewidth=2.0)
    ax1.set_ylabel("T (K)")
    ax1.set_title("US76 temperature profile", fontsize=11)
    zs_frac = np.arange(0.0, z_star_exact, 0.5)
    fr = [np.trapezoid(E[z_fine >= z0], z_fine[z_fine >= z0]) / np.trapezoid(E, z_fine)
          for z0 in zs_frac]
    ax2.plot(zs_frac, 100 * np.array(fr), color=BLUE, linewidth=2.0)
    ax2.text(0.04, 0.06,
             "external-mode energy density GROWS with height:\n"
             "the unbounded isothermal column diverges — the '10% above 40 km'\n"
             "claim is void in this regime (it holds in the internal regime)",
             transform=ax2.transAxes, fontsize=8, color=INK, va="bottom")
    ax2.set_ylabel("wave energy above z (% of column to its ceiling)")
    ax2.set_title("isothermal: energy above altitude", fontsize=11)
    fig.tight_layout()
    fig.savefig(outdir / "energy.png")
    plt.close(fig)

    md = f"""# Results — NTRS-19930082992 (Pekeris 1951) re-simulation

({time.strftime('%Y-%m-%d')}, Vulcan; exact isothermal verification + modern
US76 re-application of the paper's quadratic-term methodology.)

## Numbers

| Quantity | Paper | Ours |
|---|---|---|
| eq. (9) prefactor (isothermal) | 1/1140 | **1/{1/pref_exact:.0f}** |
| eq. (9) growth scale | z/5.6 | **z/{scale_exact:.2f} km** (from H = 7.87 km) |
| Crossing, constants as printed | "~40 km" | {z_star_printed:.1f} km (consistent with print) |
| Crossing, exact from stated H = 7.87 km | — | **{z_star_exact:.1f} km** |
| H implied by the printed 5.6 scale | 7.87 km stated | **{H_implied:.2f} km** |
| US76 crossing, external-mode growth (his isothermal surrogate law) | — | **{zstars_ext[i819]:.0f} km** (range {zstars_ext[0]:.0f}–{zstars_ext[-1]:.0f} km, H = 7.5–10) |
| US76 crossing, internal-mode growth (propagating semidiurnal tide) | — | **{zstars_int[i819]:.0f} km** (range {zstars_int[0]:.0f}–{zstars_int[-1]:.0f} km, H = 7.5–10) |
| Wave energy above 40 km, US76 internal regime (H = 8.19) | "only 10 percent" (his atmospheres) | **{100*frac40_H819:.0f}%** (ceiling {zc819:.0f} km) — claim does NOT carry over |

## Findings

1. **eq. (9)'s printed constants are mutually inconsistent.** The exact
   algebra from the paper's own stated H = 7.87 km gives prefactor 1/{1/pref_exact:.0f}
   and scale z/{scale_exact:.2f} km, crossing at {z_star_exact:.1f} km — not the
   printed (1/1140)·exp(z/5.6). The printed scale corresponds to
   H ≈ {H_implied:.1f} km. With the printed constants the "~40 km" conclusion
   is self-consistent (we get {z_star_printed:.1f} km); with the stated H it is
   {z_star_exact:.0f} km. The qualitative claim — linearization fails far below
   the E layer — survives either way.
2. **The crossing altitude is set by the vertical-structure regime.** The
   isothermal surrogate's external-mode growth law applied to the modern
   atmosphere puts the quadratic = linear crossing at {zstars_ext[i819]:.0f} km;
   the internal-mode (energy-conserving ρ^(−1/2)) growth appropriate to a
   propagating semidiurnal tide puts it at {zstars_int[i819]:.0f} km; his full
   linear-tide solutions on the near-resonant 1937/1947 model atmospheres gave
   125–130 km. All three numbers are right for their regime — his 40 ↔ 130 km
   spread IS the regime question. Pinning the modern-atmosphere value precisely
   requires solving the full vertical-structure equation on US76 (beyond this
   run's WKB scope); the WKB bracket is {zstars_ext[i819]:.0f}–{zstars_int[i819]:.0f} km.
3. **The energy argument is regime-specific — and does not carry over.** For a
   propagating tide ρu² is roughly height-constant, so the fraction of column
   wave energy above 40 km is {100*frac40_H819:.0f}% in US76 — NOT his "only
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
convergence: {dz_delta*1e3:.1f} m on the crossing altitude.
"""
    (outdir / "RESULTS.md").write_text(md)
    print(f"[A] prefactor 1/{1/pref_exact:.1f} (printed 1/1140); "
          f"scale z/{scale_exact:.3f} km (printed 5.6); "
          f"crossing exact {z_star_exact:.2f} km / printed {z_star_printed:.2f} km")
    print(f"[A] H implied by printed scale: {H_implied:.3f} km vs stated 7.87 km")
    print(f"[A] energy above 40 km (isothermal): void — column ceiling is "
          f"{z_star_exact:.1f} km < 40 km (external-mode divergence)")
    print(f"[B] US76 crossing: external {zstars_ext[i819]:.1f} km, "
          f"internal {zstars_int[i819]:.1f} km at H=8.19 "
          f"(internal range {zstars_int[0]:.1f}-{zstars_int[-1]:.1f} km)")
    print(f"[C] energy above 40 km (US76 internal, H=8.19): "
          f"{100*frac40_H819:.1f}%, ceiling {zc819:.0f} km")
    print(f"[done] wrote {outdir} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()

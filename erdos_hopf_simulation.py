#!/usr/bin/env python3
"""
Erdős Unit Distance Problem via Hopf Fibration Geometric Embedding
==================================================================

Three-stage topological transformation pipeline:
  Stage 1 │ S³ (4D sphere)       — generate Hopf fibres via golden-angle S² spiral
  Stage 2 │ ℝ³ (Clifford Torus)  — stereographic projection 4D → 3D
  Stage 3 │ ℝ² (Phyllotaxis)     — conformal flattening; each fibre → unit circle

Analytical engine:
  KD-Tree counts all unit-distance pairs ν(n) across increasing n,
  then fits power law ν(n) ≈ c·nᵅ to measure the asymptotic exponent α.

Usage:
  python3 erdos_hopf_simulation.py              # default: 300 fibres
  python3 erdos_hopf_simulation.py --fibers 500 --theta 96 --steps 16

Math reference:
  Hopf map  η: S³ → S²   η(z₁,z₂) = (2 Re(z₁z̄₂), 2 Im(z₁z̄₂), |z₁|²−|z₂|²)
  Inverse lift (section):
      z₁ = √((1+c)/2)·e^{-iφ/2},  z₂ = √((1−c)/2)·e^{+iφ/2}
  Stereo from N=(0,0,0,1):
      X = x₁/(1−x₄),  Y = x₂/(1−x₄),  Z = x₃/(1−x₄)
  Phyllotaxis:
      rₖ = √(k+0.5),  φₖ = k·φ_gold   (φ_gold ≈ 137.508°)
"""

import argparse
import sys
import warnings

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D          # noqa: F401 – registers projection
from scipy.spatial import KDTree

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────
# GLOBAL CONSTANTS
# ─────────────────────────────────────────────────────────────

# Golden angle in radians: φ_gold = 2π(1 − 1/φ²) = π(3−√5)  ≈ 137.508°
GOLDEN_ANGLE = np.pi * (3.0 - np.sqrt(5.0))

# Tolerance for |d − 1| ≤ ε  (unit-distance test)
UNIT_EPS = 1e-5

# Colour palette (dark-space theme)
COL_BG       = "#060612"
COL_PANEL    = "#0a0a1c"
COL_FIBRE_CM = "plasma"
COL_GOLD     = "#ffd700"
COL_CYAN     = "#00e5ff"
COL_GREEN    = "#00e676"
COL_ORANGE   = "#ffab40"


# ══════════════════════════════════════════════════════════════
# STAGE 1 — HOPF FIBRATION ON S³
# ══════════════════════════════════════════════════════════════

def golden_spiral_s2(n: int) -> np.ndarray:
    """
    Distribute n points on S² with uniform area density via the golden-angle spiral.

    Returns
    -------
    pts : (n, 3) array of unit vectors
    """
    k = np.arange(n, dtype=np.float64)
    # Uniform in cos θ ensures equal-area distribution
    cos_th = 1.0 - 2.0 * (k + 0.5) / n
    sin_th = np.sqrt(np.clip(1.0 - cos_th**2, 0.0, 1.0))
    phi    = GOLDEN_ANGLE * k
    return np.stack([sin_th * np.cos(phi),
                     sin_th * np.sin(phi),
                     cos_th], axis=1)                    # (n, 3)


def hopf_lift(s2_pts: np.ndarray) -> tuple:
    """
    Lift S² → S³ via the canonical section of the Hopf map.

    Given unit vector (a, b, c) ∈ S², returns (z₁, z₂) ∈ ℂ² with
    |z₁|² + |z₂|² = 1:

        z₁ = √((1+c)/2) · e^{-iφ/2},   φ = atan2(b, a)
        z₂ = √((1−c)/2) · e^{+iφ/2}

    Parameters
    ----------
    s2_pts : (n, 3)

    Returns
    -------
    z1, z2 : complex arrays of shape (n,)
    """
    a, b, c = s2_pts[:, 0], s2_pts[:, 1], s2_pts[:, 2]
    phi = np.arctan2(b, a)
    z1  = np.sqrt(np.clip((1.0 + c) / 2.0, 0.0, 1.0)) * np.exp(-0.5j * phi)
    z2  = np.sqrt(np.clip((1.0 - c) / 2.0, 0.0, 1.0)) * np.exp(+0.5j * phi)
    return z1, z2


def build_hopf_fibers(n_fibers: int, n_theta: int) -> np.ndarray:
    """
    Generate the full Hopf fibre bundle over a golden-angle S² grid.

    For each base point (z₁₀, z₂₀) ∈ S³, the Hopf fibre is:
        { (z₁₀·e^{iθ}, z₂₀·e^{iθ}) : θ ∈ [0, 2π) }

    Each fibre is a great circle on S³ whose Hopf image is a single
    point on S².

    Returns
    -------
    fibers_s3 : (n_fibers, n_theta, 4) real coordinates in ℝ⁴
        Layout: [Re z₁, Im z₁, Re z₂, Im z₂]
    """
    s2_grid        = golden_spiral_s2(n_fibers)
    z1_base, z2_base = hopf_lift(s2_grid)              # (n_fibers,)

    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    phase = np.exp(1j * theta)                         # (n_theta,)

    # Broadcast: (n_fibers, 1) * (1, n_theta) → (n_fibers, n_theta)
    z1 = z1_base[:, None] * phase[None, :]
    z2 = z2_base[:, None] * phase[None, :]

    return np.stack([z1.real, z1.imag, z2.real, z2.imag], axis=-1)


# ══════════════════════════════════════════════════════════════
# STAGE 2 — STEREOGRAPHIC PROJECTION  S³ → ℝ³
# ══════════════════════════════════════════════════════════════

def stereo_s3_to_r3(pts4: np.ndarray, pole_guard: float = 1e-9) -> np.ndarray:
    """
    Stereographic projection from S³ ⊂ ℝ⁴ to ℝ³.

    Projection from North Pole N = (0, 0, 0, 1):
        X = x₁/(1−x₄),  Y = x₂/(1−x₄),  Z = x₃/(1−x₄)

    Coordinate layout in input: (x₁, x₂, x₃, x₄) = (Re z₁, Im z₁, Re z₂, Im z₂)

    Parameters
    ----------
    pts4 : (..., 4)

    Returns
    -------
    pts3 : (..., 3)
    """
    x1, x2, x3, x4 = (pts4[..., i] for i in range(4))
    denom = np.where(np.abs(1.0 - x4) < pole_guard, pole_guard, 1.0 - x4)
    return np.stack([x1 / denom, x2 / denom, x3 / denom], axis=-1)


# ══════════════════════════════════════════════════════════════
# STAGE 3 — CONFORMAL FLATTENING → ℝ²  (Phyllotaxis lattice)
# ══════════════════════════════════════════════════════════════

def phyllotaxis_centers(n_fibers: int) -> np.ndarray:
    """
    Map each Hopf fibre to a 2D centre via the golden-angle phyllotaxis formula.

    This is the conformal image of the loxodromic strip unwrapping:
        rₖ = √(k + 0.5)   (uniform area density)
        φₖ = k · φ_gold   (golden-angle rotation prevents radial alignment)

    Returns
    -------
    centers : (n_fibers, 2)
    """
    k   = np.arange(n_fibers, dtype=np.float64)
    r   = np.sqrt(k + 0.5)
    phi = k * GOLDEN_ANGLE
    return np.stack([r * np.cos(phi), r * np.sin(phi)], axis=1)


def sample_unit_circles(centers: np.ndarray, n_theta: int) -> np.ndarray:
    """
    Discretise the unit circle around each centre.

    Parameters
    ----------
    centers : (n, 2)
    n_theta : int — points per circle

    Returns
    -------
    pts : (n * n_theta, 2)
    """
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    cx = centers[:, 0, None] + cos_t[None, :]   # (n, n_theta)
    cy = centers[:, 1, None] + sin_t[None, :]
    return np.stack([cx, cy], axis=-1).reshape(-1, 2)


# ══════════════════════════════════════════════════════════════
# ANALYTICAL ENGINE — KD-TREE UNIT DISTANCE COUNTER
# ══════════════════════════════════════════════════════════════

def circle_circle_intersections(centers: np.ndarray) -> np.ndarray:
    """
    Compute all exact pairwise intersection points of the unit circle system.

    Two unit circles with centres C_i, C_j at distance d ∈ (0, 2) intersect at:
        M = (C_i + C_j) / 2
        h = √(1 − (d/2)²)
        perp = rot90(C_j − C_i) / d
        P₁ = M + h·perp,   P₂ = M − h·perp

    Every such P is at distance exactly 1 from both C_i and C_j.

    Returns
    -------
    ipts : (m, 2)  exact intersection points (may be empty)
    """
    tree = KDTree(centers)
    ipts = []
    for i, c in enumerate(centers):
        for j in tree.query_ball_point(c, r=2.0 - 1e-9):
            if j <= i:
                continue
            diff = centers[j] - c
            d    = np.linalg.norm(diff)
            if d < 1e-9:
                continue
            h    = np.sqrt(max(0.0, 1.0 - (d * 0.5) ** 2))
            mid  = (c + centers[j]) * 0.5
            perp = np.array([-diff[1], diff[0]]) / d
            ipts.append(mid + h * perp)
            ipts.append(mid - h * perp)
    return np.array(ipts) if ipts else np.empty((0, 2))


def build_point_set(n_fibers: int, n_theta: int) -> np.ndarray:
    """
    Construct the full planar point set from the Hopf-phyllotaxis embedding.

    Combines three layers:
      • Phyllotaxis centres         — n_fibers points
      • Discretised unit circles    — n_fibers × n_theta points
        (n_theta MUST be divisible by 6 so that intra-circle distance-1
         pairs exist: 2·sin(π·(n_theta/6)/n_theta) = 2·sin(π/6) = 1 ✓)
      • Circle-circle intersections — exact unit-distance nodes (caustic layer)

    Returns
    -------
    pts : (N, 2) unique 2D points
    """
    assert n_theta % 6 == 0, f"n_theta={n_theta} must be divisible by 6"
    centers  = phyllotaxis_centers(n_fibers)
    circ_pts = sample_unit_circles(centers, n_theta)
    isect    = circle_circle_intersections(centers)

    if len(isect):
        all_pts = np.vstack([centers, circ_pts, isect])
    else:
        all_pts = np.vstack([centers, circ_pts])

    return np.unique(np.round(all_pts, 8), axis=0)


def count_unit_pairs(pts: np.ndarray, eps: float = UNIT_EPS) -> int:
    """
    Count all distinct pairs (i, j) with  1 − ε ≤ ‖pᵢ − pⱼ‖ ≤ 1 + ε.

    Strategy:
      1. KD-Tree query_pairs retrieves every pair within radius 1 + ε  (O(n·k̄))
      2. Vectorised distance filter removes those with d < 1 − ε

    Parameters
    ----------
    pts : (n, 2)
    eps : float

    Returns
    -------
    count : int
    """
    if len(pts) < 2:
        return 0
    tree = KDTree(pts)
    candidates = tree.query_pairs(r=1.0 + eps, output_type="ndarray")
    if len(candidates) == 0:
        return 0
    d = np.linalg.norm(pts[candidates[:, 0]] - pts[candidates[:, 1]], axis=1)
    return int(np.sum(d >= 1.0 - eps))


def asymptotic_scan(
    fiber_counts: list,
    n_theta: int = 30,
    eps: float = UNIT_EPS,
) -> tuple:
    """
    Sweep over increasing fibre counts; return (n_points, nu) arrays.

    Parameters
    ----------
    fiber_counts : list[int] — values of n_fibers to evaluate
    n_theta      : int — circle sample points, must be divisible by 6
                         (default 30 = 6×5 — light enough for fast scans)
    eps          : float — unit-distance tolerance

    Returns
    -------
    n_arr  : (m,) float array — total unique points at each step
    nu_arr : (m,) float array — unit-distance pair count at each step
    """
    # Silently round up to nearest multiple of 6
    if n_theta % 6 != 0:
        n_theta = ((n_theta // 6) + 1) * 6

    n_list, nu_list = [], []
    hdr = f"{'n_fibers':>10}  {'n_points':>10}  {'ν(n)':>12}"
    print(f"\n  {hdr}\n  {'─'*38}")

    for nf in fiber_counts:
        pts_u = build_point_set(nf, n_theta)
        nu    = count_unit_pairs(pts_u, eps)
        n_list.append(len(pts_u))
        nu_list.append(nu)
        print(f"  {nf:>10d}  {len(pts_u):>10,d}  {nu:>12,d}")

    return np.asarray(n_list, dtype=float), np.asarray(nu_list, dtype=float)


def fit_power_law(n: np.ndarray, nu: np.ndarray) -> tuple:
    """
    Fit ν(n) = c · nᵅ by log-linear (OLS) regression.

    Returns
    -------
    c, alpha, r_squared : float
    """
    mask = (n > 0) & (nu > 0)
    ln, lnu = np.log(n[mask]), np.log(nu[mask])
    p    = np.polyfit(ln, lnu, 1)
    pred = np.polyval(p, ln)
    ss_res = np.sum((lnu - pred) ** 2)
    ss_tot = np.sum((lnu - lnu.mean()) ** 2)
    r2 = float(1.0 - ss_res / max(ss_tot, 1e-15))
    return float(np.exp(p[1])), float(p[0]), r2


# ══════════════════════════════════════════════════════════════
# VISUALISATION
# ══════════════════════════════════════════════════════════════

def _ax_style(ax, title="", xlabel="", ylabel="",
              tc="#c8c8ff", ac="#7070aa"):
    """Apply shared dark-theme styles to a 2-D axes."""
    ax.set_facecolor(COL_PANEL)
    ax.set_title(title, color=tc, fontsize=9, pad=6, fontweight="bold")
    ax.set_xlabel(xlabel, color=ac, fontsize=8)
    ax.set_ylabel(ylabel, color=ac, fontsize=8)
    ax.tick_params(colors=ac, labelsize=6)
    for sp in ax.spines.values():
        sp.set_edgecolor("#252545")


def draw_clifford_torus(ax, fibers_3d: np.ndarray, n_show: int = 140):
    """
    Render the Clifford Torus in 3D:
      • Coloured Hopf fibres as closed curves
      • Glowing golden nodes at each fibre's geometric centre
    """
    n = min(n_show, len(fibers_3d))
    cmap = plt.get_cmap(COL_FIBRE_CM)
    ax.set_facecolor(COL_PANEL)

    centres = fibers_3d[:n].mean(axis=1)          # (n, 3)

    for i in range(n):
        f = fibers_3d[i]                           # (n_theta, 3)
        col = (*cmap(i / n)[:3], 0.50)
        f_cl = np.vstack([f, f[[0]]])              # close the loop
        ax.plot(f_cl[:, 0], f_cl[:, 1], f_cl[:, 2],
                color=col, linewidth=0.75)

    # Multi-pass glow on fibre centres
    for sz, al in [(110, 0.06), (45, 0.28), (14, 0.90)]:
        ax.scatter(*centres.T, c=COL_GOLD, s=sz, alpha=al,
                   zorder=6, depthshade=False, linewidths=0)

    ax.set_title(
        "Clifford Torus  ℝ³\n(Hopf Fibres via Stereo 4D→3D)",
        color="#b0b0ff", fontsize=9, pad=6, fontweight="bold",
    )
    for lbl, fn in zip(["X", "Y", "Z"],
                        [ax.set_xlabel, ax.set_ylabel, ax.set_zlabel]):
        fn(lbl, color="#4444aa", fontsize=7)
    ax.tick_params(colors="#33339a", labelsize=5)
    ax.xaxis.pane.fill = ax.yaxis.pane.fill = ax.zaxis.pane.fill = False
    ax.grid(False)
    ax.xaxis.pane.set_edgecolor("#1a1a3a")
    ax.yaxis.pane.set_edgecolor("#1a1a3a")
    ax.zaxis.pane.set_edgecolor("#1a1a3a")


def draw_phyllotaxis(ax, centers: np.ndarray,
                     n_theta: int = 256, n_show: int = None):
    """
    Render the 2D phyllotaxis unit-circle network:
      • Coloured unit circles (r = 1) around each phyllotaxis centre
      • Glowing caustic nodes at circle-circle intersection points
      • Centre scatter with plasma colour map
    """
    n = len(centers) if n_show is None else min(n_show, len(centers))
    ctr = centers[:n]
    cmap = plt.get_cmap("viridis")

    # ── Draw unit circles ────────────────────────────────────
    t = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    cos_t, sin_t = np.cos(t), np.sin(t)

    for i, (cx, cy) in enumerate(ctr):
        col = (*cmap(i / n)[:3], 0.28)
        ax.plot(cx + cos_t, cy + sin_t, color=col, linewidth=0.40, zorder=2)

    # ── Centre scatter (glow layers) ─────────────────────────
    c_vals = np.arange(n)
    for sz, al in [(90, 0.05), (28, 0.22), (6, 0.95)]:
        ax.scatter(ctr[:, 0], ctr[:, 1],
                   c=c_vals, cmap="plasma", s=sz, alpha=al,
                   zorder=4, linewidths=0)

    # ── Caustic nodes: circle-circle intersection points ─────
    tree = KDTree(ctr)
    caustic = []
    check_n = min(n, 350)
    for i, c in enumerate(ctr[:check_n]):
        for j in tree.query_ball_point(c, r=2.0 - 1e-6):
            if j <= i:
                continue
            d = np.linalg.norm(ctr[j] - c)
            if d < 1e-3 or d > 2.0 - 1e-6:
                continue
            # Two unit circles with centre distance d intersect at two points:
            #   midpoint M, perpendicular offset h = √(1 − (d/2)²)
            h   = np.sqrt(max(0.0, 1.0 - (d / 2.0) ** 2))
            mid = (c + ctr[j]) / 2.0
            perp = np.array([-(ctr[j][1] - c[1]), ctr[j][0] - c[0]]) / d
            caustic.extend([mid + h * perp, mid - h * perp])

    if caustic:
        cp = np.array(caustic)
        ax.scatter(cp[:, 0], cp[:, 1],
                   c=COL_GOLD, s=4, alpha=0.40, zorder=5, linewidths=0)

    r_max = np.sqrt(n) * 1.30
    ax.set_xlim(-r_max, r_max)
    ax.set_ylim(-r_max, r_max)
    ax.set_aspect("equal")
    _ax_style(
        ax,
        title="Phyllotaxis Network — Unit Circles  ℝ²\n(Conformal Flattening of Clifford Torus)",
        xlabel="x", ylabel="y",
        tc=COL_GREEN, ac="#50aa50",
    )
    ax.grid(True, color="#0c180c", linewidth=0.30)


def draw_asymptote(ax, n_arr: np.ndarray, nu_arr: np.ndarray,
                   c: float, alpha: float, r2: float):
    """
    Log-log plot of ν(n) with:
      • fitted power law  ν ≈ c·nᵅ
      • Szemerédi-Trotter upper bound  O(n^{4/3})
      • Erdős square-grid baseline  n^{1 + c/ln ln n}
    """
    if len(n_arr) < 2:
        return

    n_sm = np.logspace(np.log10(n_arr.min()), np.log10(n_arr.max()), 400)

    # Measured
    ax.scatter(n_arr, nu_arr, color=COL_CYAN, s=40, zorder=6,
               label="Measured ν(n)")
    ax.plot(n_arr, nu_arr, color=COL_CYAN, alpha=0.35, linewidth=0.9)

    # Power-law fit
    ax.plot(n_sm, c * n_sm ** alpha,
            "--", color=COL_ORANGE, linewidth=2.2, zorder=7,
            label=rf"Fit: $c\,n^{{{alpha:.4f}}}$  ($R^2={r2:.3f}$)")

    # Szemerédi-Trotter ceiling
    ax.plot(n_sm, 0.5 * n_sm ** (4.0 / 3.0),
            "-.", color="#ff4040", linewidth=1.2, alpha=0.75,
            label=r"S-T bound:  $O(n^{4/3})$")

    # Erdős square-grid baseline  (exponent → 1 as n → ∞)
    with np.errstate(divide="ignore", invalid="ignore"):
        log_log_n = np.log(np.log(np.maximum(n_sm, np.e**np.e)))
        erdos_exp = 1.0 + 1.0 / log_log_n
    ax.plot(n_sm, n_sm ** erdos_exp,
            ":", color="#ff8030", linewidth=1.3, alpha=0.75,
            label=r"Grid: $n^{1+c/\ln\!\ln n}$")

    ax.set_xscale("log")
    ax.set_yscale("log")
    _ax_style(
        ax,
        title=f"Asymptotic Growth  α ≈ {alpha:.4f}",
        xlabel="n  (total points)",
        ylabel="ν(n)  (unit-distance pairs)",
        tc="#ffff88", ac="#aaaa55",
    )
    ax.legend(fontsize=6.5, facecolor="#101025", edgecolor="#3a3a6a",
              labelcolor="white", framealpha=0.88)
    ax.grid(True, which="both", color="#18180a", linewidth=0.30)


def draw_dist_histogram(ax, pts: np.ndarray, sample_n: int = 700):
    """
    Histogram of sampled pairwise distances highlighting the d = 1 peak.
    """
    rng = np.random.default_rng(0)
    idx = rng.choice(len(pts), min(sample_n, len(pts)), replace=False)
    smp = pts[idx]

    tree = KDTree(smp)
    dists = []
    for i, p in enumerate(smp):
        for j in tree.query_ball_point(p, r=2.6):
            if j > i:
                dists.append(np.linalg.norm(smp[j] - p))

    if not dists:
        ax.text(0.5, 0.5, "No distances", ha="center", va="center",
                transform=ax.transAxes, color="white")
        return

    dists = np.asarray(dists)
    bins = np.linspace(0.0, 2.6, 110)
    ax.hist(dists, bins=bins, color="#2a5a9a", edgecolor="none", alpha=0.85)
    ax.axvspan(0.90, 1.10, alpha=0.22, color=COL_GOLD)
    ax.axvline(1.0, color=COL_GOLD, linewidth=2.0, linestyle="--",
               label="d = 1  (unit)")
    _ax_style(
        ax,
        title="Pairwise Distance Distribution\n(unit-distance peak at d = 1)",
        xlabel="d", ylabel="count",
        tc="#88ccff", ac="#5090cc",
    )
    ax.legend(fontsize=7, facecolor="#101025", edgecolor="#3a3a6a",
              labelcolor="white", framealpha=0.88)
    ax.grid(True, color="#0a0a1a", linewidth=0.30)


# ══════════════════════════════════════════════════════════════
# MAIN SIMULATION
# ══════════════════════════════════════════════════════════════

def run_simulation(
    n_fibers: int      = 300,
    n_theta:  int      = 96,    # must be divisible by 6 (96 = 6×16)
    n_scan_steps: int  = 14,
    output: str        = "erdos_hopf_simulation.png",
    show: bool         = True,
):
    """
    Execute the full Hopf-Fibration unit-distance simulation.

    Parameters
    ----------
    n_fibers     : number of Hopf fibres (= unit circles placed in ℝ²)
    n_theta      : discretisation points per fibre / circle
                   MUST be divisible by 6 (guarantees d=1 intra-circle pairs)
    n_scan_steps : resolution of the asymptotic sweep
    output       : filename for saved PNG figure
    show         : whether to call plt.show() at the end
    """
    # Enforce divisibility-by-6 constraint
    if n_theta % 6 != 0:
        n_theta = ((n_theta // 6) + 1) * 6
        print(f"[info]  n_theta rounded up to {n_theta} (nearest 6-multiple)")

    SEP = "═" * 68
    print(f"\n{SEP}")
    print("  ERDŐS UNIT DISTANCE — HOPF FIBRATION SIMULATION")
    print(f"  n_fibers={n_fibers}  n_theta={n_theta}")
    print(SEP)

    # ── STAGE 1: Hopf fibres on S³ ───────────────────────────
    print(f"\n[1/4]  Generating {n_fibers} Hopf fibres on S³ …")
    fibers_s3 = build_hopf_fibers(n_fibers, n_theta)
    print(f"       Tensor shape: {fibers_s3.shape}  (fibres × θ-samples × ℝ⁴)")

    # ── STAGE 2: Stereographic projection S³ → ℝ³ ───────────
    print("\n[2/4]  Stereographic projection S³ → ℝ³ (Clifford Torus) …")
    fibers_3d = stereo_s3_to_r3(fibers_s3)
    fibers_3d = np.clip(fibers_3d, -25.0, 25.0)    # suppress values near North Pole
    flat3 = fibers_3d.reshape(-1, 3)
    print(f"       3D bounding box: [{flat3.min():.3f}, {flat3.max():.3f}]")

    # ── STAGE 3: Phyllotaxis lattice ─────────────────────────
    print("\n[3/4]  Conformal flattening → 2D phyllotaxis lattice …")
    centers_2d = phyllotaxis_centers(n_fibers)
    pts_2d     = build_point_set(n_fibers, n_theta)
    print(f"       Unique 2D points: {len(pts_2d):,}")
    print(f"       (circles={n_fibers}×{n_theta} + centres + ∩-points)")

    # Quick unit-distance count (capped for speed)
    cap    = min(len(pts_2d), 8000)
    nu_main = count_unit_pairs(pts_2d[:cap])
    print(f"       Unit-distance pairs (first {cap:,} pts): {nu_main:,}")

    # ── ASYMPTOTIC SCAN ───────────────────────────────────────
    print("\n[4/4]  Asymptotic scan …")
    hi = min(n_fibers, 120)
    fiber_sweep = np.unique(
        np.round(np.logspace(np.log10(8), np.log10(hi), n_scan_steps)).astype(int)
    ).tolist()
    n_arr, nu_arr = asymptotic_scan(fiber_sweep, n_theta=30)   # 30 = 6×5

    c_fit = alpha_fit = r2_fit = float("nan")
    if len(n_arr) >= 3 and np.all(nu_arr > 0):
        c_fit, alpha_fit, r2_fit = fit_power_law(n_arr, nu_arr)
        print(f"\n  Power-law fit :  ν(n) ≈ {c_fit:.5f} · n^{alpha_fit:.6f}")
        print(f"  R²             = {r2_fit:.5f}")
        print(f"  S-T upper bound: 4/3 ≈ {4/3:.6f}")
        print(f"  Erdős baseline : α → 1.0000 as n → ∞")
        margin = alpha_fit - 1.0
        print(f"  Exponent gain  : Δα = {margin:+.6f} above linear")

    # ── BUILD FIGURE ──────────────────────────────────────────
    print("\n[Render]  Building split-screen figure …")
    fig = plt.figure(figsize=(22, 12), facecolor=COL_BG)
    fig.suptitle(
        "Erdős Unit Distance Problem — Hopf Fibration Geometric Embedding\n"
        r"$S^3\;[\text{Hopf}]\;\longrightarrow\;"
        r"\mathbb{R}^3\;[\text{Stereo}]\;\longrightarrow\;"
        r"\mathbb{R}^2\;[\text{Conformal}]$",
        color="#d8d8ff", fontsize=12, fontweight="bold", y=0.995,
    )

    gs = gridspec.GridSpec(
        2, 3, figure=fig,
        left=0.03, right=0.97, top=0.91, bottom=0.05,
        hspace=0.40, wspace=0.28,
    )
    ax_3d   = fig.add_subplot(gs[:, 0], projection="3d")
    ax_2d   = fig.add_subplot(gs[:, 1])
    ax_asym = fig.add_subplot(gs[0, 2])
    ax_hist = fig.add_subplot(gs[1, 2])

    draw_clifford_torus(ax_3d, fibers_3d, n_show=min(n_fibers, 150))
    draw_phyllotaxis(ax_2d, centers_2d, n_theta=256,
                     n_show=min(n_fibers, 500))
    if not np.isnan(alpha_fit):
        draw_asymptote(ax_asym, n_arr, nu_arr, c_fit, alpha_fit, r2_fit)
    draw_dist_histogram(ax_hist, pts_2d[:min(len(pts_2d), 6000)])

    # ── Summary annotation ────────────────────────────────────
    exp_str = f"{alpha_fit:.5f}" if not np.isnan(alpha_fit) else "N/A"
    r2_str  = f"{r2_fit:.4f}"   if not np.isnan(r2_fit)    else "N/A"
    ann = (
        f"n_fibers  = {n_fibers:>6d}\n"
        f"n_points  = {len(pts_2d):>8,d}\n"
        f"ν(n) est  = {nu_main:>8,d}\n"
        f"──────────────────\n"
        f"α (fit)   = {exp_str:>8s}\n"
        f"R²        = {r2_str:>8s}\n"
        f"S-T bound = {'1.3333':>8s}\n"
        f"Grid base = {'1.0000':>8s}"
    )
    fig.text(
        0.977, 0.03, ann, ha="right", va="bottom",
        color="#88ff88", fontsize=7.5, fontfamily="monospace",
        bbox=dict(boxstyle="round,pad=0.6", facecolor="#070f07",
                  edgecolor="#205520", alpha=0.92),
    )

    plt.savefig(output, dpi=150, bbox_inches="tight",
                facecolor=COL_BG, edgecolor="none")
    print(f"\n[Done]  Figure saved → {output}")

    if show:
        plt.show()

    return fig, n_arr, nu_arr, alpha_fit


# ══════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Erdős Unit Distance — Hopf Fibration Simulation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--fibers",  type=int, default=300,
                   help="Number of Hopf fibres")
    p.add_argument("--theta",   type=int, default=96,
                   help="Sample points per fibre / unit circle (must be divisible by 6)")
    p.add_argument("--steps",   type=int, default=14,
                   help="Steps in asymptotic sweep")
    p.add_argument("--output",  type=str, default="erdos_hopf_simulation.png",
                   help="Output PNG filename")
    p.add_argument("--no-show", action="store_true",
                   help="Skip plt.show() (useful in headless environments)")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    run_simulation(
        n_fibers     = args.fibers,
        n_theta      = args.theta,
        n_scan_steps = args.steps,
        output       = args.output,
        show         = not args.no_show,
    )

"""
Stream-function BEM solver for the hermetic Meissner-sheath toroidal axion
pickup.  Implements the unit-current boundary-value problem laid out in
cc_task_spec.md and Sec. 2-3 of toroid_optimization.pdf.

Conventions
-----------
* Length in units of major radius R = 1, mu_0 = 1, unit current I = 1.
* Sheath bounds V_mag = {1 < s < 1+beta, 0 < z < alpha, 0 < phi < 2*pi}.
* B_0 = (R/s) phi_hat inside V_mag (we set R = B_max = 1; physical
  prefactors are restored in the crossover-R_t calculation).
* "Outward from V_mag" normals:
      A (inner cyl, s=1):    n = -s_hat
      B (top cap,  z=alpha): n = +z_hat
      C (outer cyl, s=1+b):  n = +s_hat
      D (bottom cap, z=0):   n = -z_hat
* Sheath = topological torus.  Slit ribbon at |phi| < w/2 around the
  meridional rectangle removes a topological annulus, leaving a
  topological cylinder.  Stream function psi is single-valued on this
  cylinder, with Dirichlet data
      psi = +1/2  on phi = +w/2  (right edge),
      psi = -1/2  on phi = -w/2  (left edge),
  i.e. unit current driven between the two slit edges.
* External return wire is parameterised as a polyline; carries +1 A
  flowing from terminal_+ to terminal_-.

Energy formulation
------------------
The induced surface current K (subject to Meissner BC) minimises the
total magnetic energy
    W = 1/2 int K . A_total dS
       = (mu_0/8pi) [ int int K(r).K(r')/|r-r'| dS dS'
                    + 2 int int K(r) . dl'(r')/|r-r'| dS dl'_y
                    + L_wire_self ]
subject to the Dirichlet data on psi.  The variational equation
    delta W / delta psi_i = 0  for interior i
is equivalent to n_hat . B = 0 collocation.

Implementation notes
--------------------
The sheath is meshed by a structured (chi, t) grid with chi = phi - w/2
running in [0, 2*pi - w] and t = arc-length along the meridional
rectangle A -> B -> C -> D -> A.  Rectangle elements; bilinear nodal psi
with centroid-evaluated piecewise-constant surface current per element.
Mutual inductance kernel is evaluated via centroid quadrature with
analytic self-panel correction (flat-strip approximation).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from dataclasses import dataclass, field, asdict
from typing import Tuple, List, Dict, Optional

import numpy as np
from scipy.linalg import solve as la_solve
from scipy.linalg import lu_factor, lu_solve


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def cross_section_param(t: np.ndarray, alpha: float, beta: float):
    """Map cross-section arc-length t in [0, 2*alpha+2*beta] to (s, z, seg).

    Segments A,B,C,D have ids 0,1,2,3.
    """
    t = np.asarray(t)
    s = np.empty_like(t, dtype=float)
    z = np.empty_like(t, dtype=float)
    seg = np.empty_like(t, dtype=int)

    mask_a = t <= alpha
    mask_b = (t > alpha) & (t <= alpha + beta)
    mask_c = (t > alpha + beta) & (t <= 2 * alpha + beta)
    mask_d = t > 2 * alpha + beta

    s[mask_a] = 1.0
    z[mask_a] = t[mask_a]
    seg[mask_a] = 0

    s[mask_b] = 1.0 + (t[mask_b] - alpha)
    z[mask_b] = alpha
    seg[mask_b] = 1

    s[mask_c] = 1.0 + beta
    z[mask_c] = alpha - (t[mask_c] - alpha - beta)
    seg[mask_c] = 2

    s[mask_d] = 1.0 + beta - (t[mask_d] - 2 * alpha - beta)
    z[mask_d] = 0.0
    seg[mask_d] = 3

    return s, z, seg


def segment_normal(seg: np.ndarray) -> np.ndarray:
    """Outward (from V_mag) normal in cylindrical components (n_s, n_phi, n_z)."""
    n = np.zeros((seg.size, 3))
    n[seg == 0] = [-1.0, 0.0, 0.0]   # A: n = -s_hat
    n[seg == 1] = [0.0, 0.0, 1.0]    # B: n = +z_hat
    n[seg == 2] = [1.0, 0.0, 0.0]    # C: n = +s_hat
    n[seg == 3] = [0.0, 0.0, -1.0]   # D: n = -z_hat
    return n


def cyl_to_cart(s: np.ndarray, phi: np.ndarray, z: np.ndarray):
    """Convert cylindrical (s, phi, z) -> Cartesian (x, y, z)."""
    return s * np.cos(phi), s * np.sin(phi), z


def cyl_vec_to_cart(vs, vphi, vz, phi):
    """Convert (s,phi,z)-components of a vector at azimuth phi to Cartesian."""
    c, s_ = np.cos(phi), np.sin(phi)
    vx = vs * c - vphi * s_
    vy = vs * s_ + vphi * c
    return vx, vy, vz


# ---------------------------------------------------------------------------
# Mesh
# ---------------------------------------------------------------------------

@dataclass
class Mesh:
    alpha: float
    beta: float
    w: float
    chi_nodes: np.ndarray          # (n_chi+1,)
    t_nodes: np.ndarray            # (n_t+1,)
    chi_centers: np.ndarray        # (n_chi,)
    t_centers: np.ndarray          # (n_t,)
    n_chi: int
    n_t: int
    # Per-element fields, shape (n_chi*n_t,):
    elem_s: np.ndarray
    elem_phi: np.ndarray
    elem_z: np.ndarray
    elem_seg: np.ndarray
    elem_area: np.ndarray
    elem_normal_cart: np.ndarray   # (n_el, 3)
    elem_xyz: np.ndarray           # (n_el, 3)  Cartesian centroid
    elem_dchi: np.ndarray          # (n_el,)
    elem_dt: np.ndarray            # (n_el,)
    # Node indexing helpers
    # Total node count = (n_chi+1)*(n_t_unique). t is periodic around the
    # meridional rectangle, so node row t_node = n_t corresponds to t_node=0.
    n_t_unique: int                # = n_t (periodic in t)
    n_chi_nodes: int               # = n_chi + 1
    n_nodes: int
    # Dirichlet info: indices of left-edge (chi=0) and right-edge (chi=Phi)
    left_node_idx: np.ndarray
    right_node_idx: np.ndarray
    interior_node_idx: np.ndarray
    # Quad element corner node indices, shape (n_el, 4) in order
    # (chi_lo,t_lo), (chi_hi,t_lo), (chi_hi,t_hi), (chi_lo,t_hi)
    elem_node_idx: np.ndarray
    # Per-element basis: contribution of each corner node's psi to K_t and K_chi
    # at the element centroid.  K_t and K_chi are in the local (chi,t) frame.
    # Provided as factors d_psi/d_chi and d_psi/d_t at centroid as a linear
    # function of the 4 corner psi values.
    grad_chi_coef: np.ndarray      # (n_el, 4) -> d psi/d chi at centroid
    grad_t_coef: np.ndarray        # (n_el, 4) -> d psi/d t at centroid
    # Local frame -> global cylindrical conversion for K
    # Final K (at centroid) in cyl components (K_s, K_phi, K_z) is a linear
    # function of d psi/d chi and d psi/d t.  We store coefficients
    # K_cyl = M_chi * (d psi/d chi) + M_t * (d psi/d t)
    K_chi_coef_cyl: np.ndarray     # (n_el, 3)
    K_t_coef_cyl: np.ndarray       # (n_el, 3)
    K_chi_coef_cart: np.ndarray    # (n_el, 3)
    K_t_coef_cart: np.ndarray      # (n_el, 3)


def adaptive_chi_nodes(w: float,
                       fine_factor: float = 1.0 / 3.0,
                       n_total_hint: int = 80) -> np.ndarray:
    """Adaptive azimuthal node positions on chi in [0, 2*pi - w].

    Refinement bands (spec):
      |dphi from slit| < 3w/2:   element size <= w/3
      3w/2 < |dphi|  < 6w/2:    element size <= w
      else:                       element size <= 0.1 rad

    chi = 0 is the +w/2 slit edge; chi = 2pi - w is the -w/2 slit edge.
    Distance from slit (in phi) for a chi value is min(chi, 2*pi - w - chi).
    Build nodes piecewise.
    """
    big = 2 * np.pi - w
    far_h = 0.1
    band_inner = 3 * w / 2.0
    band_outer = 6 * w / 2.0  # = 3w
    h_inner = w * fine_factor
    h_mid = w
    h_far = far_h

    # Build node list as union of stepped regions, both ends.
    nodes = [0.0]

    def add_block(start, end, h):
        # Number of intervals
        L = end - start
        if L <= 0:
            return
        n = max(1, int(np.ceil(L / h)))
        for k in range(1, n + 1):
            nodes.append(start + k * L / n)

    # Distance from the left edge.
    # Sequence: [0, band_inner] -> h_inner
    #           [band_inner, band_outer] -> h_mid
    #           [band_outer, big - band_outer] -> h_far
    #           [big - band_outer, big - band_inner] -> h_mid
    #           [big - band_inner, big] -> h_inner
    add_block(0.0, min(band_inner, big), h_inner)
    add_block(min(band_inner, big), min(band_outer, big), h_mid)
    far_lo = min(band_outer, big)
    far_hi = max(big - band_outer, far_lo)
    add_block(far_lo, far_hi, h_far)
    far_hi2 = max(big - band_inner, far_hi)
    add_block(far_hi, far_hi2, h_mid)
    add_block(far_hi2, big, h_inner)

    nodes = np.unique(np.array(nodes))
    # Force endpoint exactness
    nodes[0] = 0.0
    nodes[-1] = big
    return nodes


def build_mesh(alpha: float,
               beta: float,
               w: float = 0.03,
               n_t_per_seg: int = 20,
               chi_nodes_override: Optional[np.ndarray] = None) -> Mesh:
    """Build the structured (chi, t) sheath mesh."""
    if chi_nodes_override is not None:
        chi_nodes = np.asarray(chi_nodes_override, dtype=float)
    else:
        chi_nodes = adaptive_chi_nodes(w)
    n_chi = chi_nodes.size - 1

    # t mesh: uniform within each of four segments, joined at corners.
    t_a = np.linspace(0.0, alpha, n_t_per_seg + 1)
    t_b = np.linspace(alpha, alpha + beta, n_t_per_seg + 1)
    t_c = np.linspace(alpha + beta, 2 * alpha + beta, n_t_per_seg + 1)
    t_d = np.linspace(2 * alpha + beta, 2 * alpha + 2 * beta, n_t_per_seg + 1)
    t_nodes = np.unique(np.concatenate([t_a, t_b, t_c, t_d]))
    n_t = t_nodes.size - 1  # 4 * n_t_per_seg

    chi_centers = 0.5 * (chi_nodes[:-1] + chi_nodes[1:])
    t_centers = 0.5 * (t_nodes[:-1] + t_nodes[1:])

    # Per-element fields: row-major (chi, t).
    CHI, T = np.meshgrid(chi_centers, t_centers, indexing="ij")
    dchi = np.diff(chi_nodes)
    dt = np.diff(t_nodes)
    DCHI, DT = np.meshgrid(dchi, dt, indexing="ij")

    elem_chi = CHI.flatten()
    elem_t = T.flatten()
    elem_dchi = DCHI.flatten()
    elem_dt = DT.flatten()
    elem_phi = elem_chi + w / 2.0  # actual azimuth

    s_cen, z_cen, seg_cen = cross_section_param(elem_t, alpha, beta)
    n_cyl = segment_normal(seg_cen)
    # Element area: depends on segment.  Cyl segments: dA = s dphi dz = s_cen * dchi * dt_z
    # Top/bottom caps: dA = s ds dphi = s_cen * dchi * dt_s.  Either way dA = s_cen * dchi * dt.
    elem_area = s_cen * elem_dchi * elem_dt

    x_cen, y_cen, _ = cyl_to_cart(s_cen, elem_phi, z_cen)
    xyz = np.stack([x_cen, y_cen, z_cen], axis=1)

    # Cartesian normal: rotate (n_s, 0, n_z) by phi.  n_phi = 0 always here.
    nx, ny, nz = cyl_vec_to_cart(n_cyl[:, 0], n_cyl[:, 1], n_cyl[:, 2], elem_phi)
    n_cart = np.stack([nx, ny, nz], axis=1)

    # Node indexing.  Topology: chi axis open (Dirichlet on ends),
    # t axis periodic (so node row t_node = n_t identified with 0).
    n_t_unique = n_t
    n_chi_nodes = n_chi + 1
    n_nodes = n_chi_nodes * n_t_unique

    def node_id(ic, it):
        return ic * n_t_unique + (it % n_t_unique)

    # Boundary node indices
    left_idx = np.array([node_id(0, it) for it in range(n_t_unique)], dtype=int)
    right_idx = np.array([node_id(n_chi, it) for it in range(n_t_unique)], dtype=int)
    bnd_idx = np.unique(np.concatenate([left_idx, right_idx]))
    all_idx = np.arange(n_nodes)
    interior_idx = np.setdiff1d(all_idx, bnd_idx)

    # Element corner node indices.  Order: (ic,it), (ic+1,it), (ic+1,it+1), (ic,it+1).
    elem_corners = np.zeros((n_chi * n_t, 4), dtype=int)
    for ic in range(n_chi):
        for it in range(n_t):
            e = ic * n_t + it
            elem_corners[e, 0] = node_id(ic, it)
            elem_corners[e, 1] = node_id(ic + 1, it)
            elem_corners[e, 2] = node_id(ic + 1, it + 1)
            elem_corners[e, 3] = node_id(ic, it + 1)

    # Bilinear ψ derivatives at element centroid:
    # ψ(ξ,η) = (1-ξ)(1-η) p0 + ξ(1-η) p1 + ξη p2 + (1-ξ)η p3
    # with ξ = (chi - chi_lo)/dchi, η = (t - t_lo)/dt
    # At centroid ξ = η = 1/2:
    #   ∂ψ/∂chi = ( -p0 + p1 + p2 - p3 ) / (2 dchi)
    #   ∂ψ/∂t   = ( -p0 - p1 + p2 + p3 ) / (2 dt)
    grad_chi_coef = np.zeros((n_chi * n_t, 4))
    grad_t_coef = np.zeros((n_chi * n_t, 4))
    grad_chi_coef[:, 0] = -1.0 / (2.0 * elem_dchi)
    grad_chi_coef[:, 1] = +1.0 / (2.0 * elem_dchi)
    grad_chi_coef[:, 2] = +1.0 / (2.0 * elem_dchi)
    grad_chi_coef[:, 3] = -1.0 / (2.0 * elem_dchi)
    grad_t_coef[:, 0] = -1.0 / (2.0 * elem_dt)
    grad_t_coef[:, 1] = -1.0 / (2.0 * elem_dt)
    grad_t_coef[:, 2] = +1.0 / (2.0 * elem_dt)
    grad_t_coef[:, 3] = +1.0 / (2.0 * elem_dt)

    # Stream-function -> K (cyl components), with global convention
    #     K = ∇_S ψ × n̂.
    # The spec's formulas K_phi = -∂_z ψ, K_z = (1/s) ∂_phi ψ on "a cylinder"
    # implicitly assume n̂ = -ŝ (inner cylinder).  On the outer cylinder
    # (n̂ = +ŝ) the meridional sign flips.  Using K = ∇_S ψ × n̂ everywhere
    # gives a single ψ that is *globally continuous* across the four
    # segments and yields a current loop that circulates poloidally
    # around the cross-section.
    #
    # With chi = phi - w/2, ∂_chi = ∂_phi, and t increasing A->B->C->D:
    #
    # Seg A (inner cyl, n̂ = -ŝ, t̂ = +ẑ):
    #     K_z   = +(1/s) ∂_chi ψ,   K_phi = -∂_t ψ
    # Seg B (top cap,   n̂ = +ẑ, t̂ = +ŝ):
    #     K_s   = +(1/s) ∂_chi ψ,   K_phi = -∂_t ψ
    # Seg C (outer cyl, n̂ = +ŝ, t̂ = -ẑ):
    #     K_z   = -(1/s) ∂_chi ψ,   K_phi = -∂_t ψ
    # Seg D (bottom cap,n̂ = -ẑ, t̂ = -ŝ):
    #     K_s   = -(1/s) ∂_chi ψ,   K_phi = -∂_t ψ
    #
    # Note the meridional component K . t̂ = (1/s) ∂_chi ψ on *all* four
    # segments, so current is conserved across corner seams.  The
    # toroidal component K_phi = -∂_t ψ is also globally consistent.

    K_chi_coef_cyl = np.zeros((n_chi * n_t, 3))   # contribution from ∂_chi ψ
    K_t_coef_cyl = np.zeros((n_chi * n_t, 3))     # contribution from ∂_t ψ

    for seg_id, idx in [(s, np.where(seg_cen == s)[0]) for s in range(4)]:
        if idx.size == 0:
            continue
        s_here = s_cen[idx]
        if seg_id == 0:        # inner cyl
            K_chi_coef_cyl[idx, 2] = +1.0 / s_here   # K_z = +(1/s) d_chi
            K_t_coef_cyl[idx, 1] = -1.0              # K_phi = -d_t
        elif seg_id == 1:      # top cap
            K_chi_coef_cyl[idx, 0] = +1.0 / s_here   # K_s = +(1/s) d_chi
            K_t_coef_cyl[idx, 1] = -1.0              # K_phi = -d_t
        elif seg_id == 2:      # outer cyl
            K_chi_coef_cyl[idx, 2] = -1.0 / s_here   # K_z = -(1/s) d_chi
            K_t_coef_cyl[idx, 1] = -1.0              # K_phi = -d_t
        elif seg_id == 3:      # bottom cap
            K_chi_coef_cyl[idx, 0] = -1.0 / s_here   # K_s = -(1/s) d_chi
            K_t_coef_cyl[idx, 1] = -1.0              # K_phi = -d_t

    # Convert these to Cartesian per element.
    cos_phi = np.cos(elem_phi)
    sin_phi = np.sin(elem_phi)

    def cyl3_to_cart3(vec_cyl, phi_):
        # vec_cyl: (n, 3) in (s, phi, z) basis at element azimuth phi_.
        vs, vphi, vz = vec_cyl[:, 0], vec_cyl[:, 1], vec_cyl[:, 2]
        c, s_ = np.cos(phi_), np.sin(phi_)
        vx = vs * c - vphi * s_
        vy = vs * s_ + vphi * c
        return np.stack([vx, vy, vz], axis=1)

    K_chi_coef_cart = cyl3_to_cart3(K_chi_coef_cyl, elem_phi)
    K_t_coef_cart = cyl3_to_cart3(K_t_coef_cyl, elem_phi)

    return Mesh(
        alpha=alpha, beta=beta, w=w,
        chi_nodes=chi_nodes, t_nodes=t_nodes,
        chi_centers=chi_centers, t_centers=t_centers,
        n_chi=n_chi, n_t=n_t,
        elem_s=s_cen, elem_phi=elem_phi, elem_z=z_cen,
        elem_seg=seg_cen, elem_area=elem_area,
        elem_normal_cart=n_cart, elem_xyz=xyz,
        elem_dchi=elem_dchi, elem_dt=elem_dt,
        n_t_unique=n_t_unique, n_chi_nodes=n_chi_nodes, n_nodes=n_nodes,
        left_node_idx=left_idx, right_node_idx=right_idx,
        interior_node_idx=interior_idx,
        elem_node_idx=elem_corners,
        grad_chi_coef=grad_chi_coef, grad_t_coef=grad_t_coef,
        K_chi_coef_cyl=K_chi_coef_cyl, K_t_coef_cyl=K_t_coef_cyl,
        K_chi_coef_cart=K_chi_coef_cart, K_t_coef_cart=K_t_coef_cart,
    )


# ---------------------------------------------------------------------------
# K from nodal psi
# ---------------------------------------------------------------------------

def K_from_psi(psi: np.ndarray, mesh: Mesh) -> np.ndarray:
    """Centroid Cartesian K (n_el, 3) given nodal psi (n_nodes,)."""
    # ∂ψ/∂chi and ∂ψ/∂t at centroids
    corners = psi[mesh.elem_node_idx]                    # (n_el, 4)
    dpsi_dchi = np.einsum("ij,ij->i", mesh.grad_chi_coef, corners)
    dpsi_dt = np.einsum("ij,ij->i", mesh.grad_t_coef, corners)
    K = (mesh.K_chi_coef_cart * dpsi_dchi[:, None]
         + mesh.K_t_coef_cart * dpsi_dt[:, None])
    return K


# ---------------------------------------------------------------------------
# External return wire (polyline filament)
# ---------------------------------------------------------------------------

@dataclass
class WirePath:
    points: np.ndarray   # (n_pt, 3) Cartesian
    radius: float = 0.01
    current: float = 1.0

    @property
    def n_seg(self) -> int:
        return self.points.shape[0] - 1


def empty_wire(r_wire: float = 0.01) -> WirePath:
    """Trivial wire with no segments: produces zero A_wire and B_wire
    everywhere, and zero L_wire_self.  Used for the sheath-only patch
    described in wire_removal_patch.md.

    With this wire, the linear system has:
      - Slit Dirichlet data (psi = +/- 1/2 on the two slit edges) as the
        only source driving non-trivial current,
      - No wire-coupling cross term in the RHS,
      - L_p = sheath self-energy only,
      - V_eff = sheath A_phi integrated over V_mag (no wire contribution).
    """
    return WirePath(points=np.zeros((1, 3)), radius=r_wire)


def production_wire(alpha: float, w: float = 0.03,
                    r_wire: float = 0.01) -> WirePath:
    """Wire path from spec.  5 straight Cartesian segments:
      P1 = (1, +w/2, alpha/2)  terminal +
      P2 = (10, +w/2, alpha/2) after radial out
      P3 = (10, +w/2, -2)      after axial down at s=10
      P4 = (0, 0, -2)          after radial in at z=-2 to axis
      P5 = (0, 0, alpha/2)     after axial up along z-axis
      P6 = (1, -w/2, alpha/2)  terminal - (final radial out)
    Note P3 -> P4 and P5 -> P6 are straight Cartesian lines that cross
    from one phi-side to the central axis; physically the wire has finite
    radius r_wire (regularizer for log singularities)."""
    plus = +w / 2.0
    minus = -w / 2.0
    z_top = alpha / 2.0
    pts_cyl = [
        (1.0, plus, z_top),
        (10.0, plus, z_top),
        (10.0, plus, -2.0),
        (0.0, 0.0, -2.0),
        (0.0, 0.0, z_top),
        (1.0, minus, z_top),
    ]
    cart = []
    for s_, p_, z_ in pts_cyl:
        cart.append([s_ * np.cos(p_), s_ * np.sin(p_), z_])
    return WirePath(points=np.array(cart, dtype=float), radius=r_wire)


def linked_wire(alpha: float, w: float = 0.03,
                r_wire: float = 0.01) -> WirePath:
    """Wire that links the major-axis cycle once -- the readout wire
    threads the central hole.  Used in XC1.  Goes from + terminal, out
    radially, up over the top, through the central hole at s = 0.3
    (passes through the magnet bore midplane), down to underneath, then
    back out radially.  The key topological feature is that the return
    path threads the central hole exactly once."""
    plus = +w / 2.0
    minus = -w / 2.0
    z_top_term = alpha / 2.0
    s_hole = 0.3
    pts_cyl = [
        (1.0, plus, z_top_term),
        (10.0, plus, z_top_term),
        (10.0, plus, alpha + 2.0),
        (s_hole, plus, alpha + 2.0),
        (s_hole, plus, -2.0),
        (10.0, plus, -2.0),
        (10.0, minus, -2.0),
        (10.0, minus, z_top_term),
        (1.0, minus, z_top_term),
    ]
    cart = []
    for s_, p_, z_ in pts_cyl:
        cart.append([s_ * np.cos(p_), s_ * np.sin(p_), z_])
    return WirePath(points=np.array(cart, dtype=float), radius=r_wire)


def wire_segments(wire: WirePath):
    """Return (starts, ends, dl) Cartesian arrays for each segment."""
    starts = wire.points[:-1]
    ends = wire.points[1:]
    return starts, ends, ends - starts


def B_from_wire(rs: np.ndarray, wire: WirePath) -> np.ndarray:
    """Biot-Savart B from finite straight segments.

    For a straight segment from a to b carrying current I in direction
    (b-a)/|b-a|, the B-field at point r is
        B = (mu_0 I / 4 pi) * (cos theta1 - cos theta2)/d  * t_hat x r_hat
    where d is the perpendicular distance, theta1, theta2 are angles at
    endpoints, t_hat is segment direction, r_hat is closest-approach
    direction from the line to r.

    We use the closed-form vector formula
        B = (mu_0 I / 4 pi) * (t_hat x r_a)/|r_a|^2 ... etc.
    via the standard finite-wire result:
        B = (mu_0 I / 4 pi |R|) * (r_b/|r_b| - r_a/|r_a|) . (t_hat) ... no.

    We just use the well-known finite-segment formula:
        Let r_a = r - a, r_b = r - b.
        Then B = (mu_0 I / 4 pi) * (r_a x r_b) /
                 ( |r_a| |r_b| (|r_a||r_b| + r_a . r_b) ).
    """
    rs = np.atleast_2d(rs)               # (n, 3)
    starts, ends, _ = wire_segments(wire)
    B = np.zeros_like(rs)
    for a, b in zip(starts, ends):
        ra = rs - a[None, :]
        rb = rs - b[None, :]
        ra_n = np.linalg.norm(ra, axis=1)
        rb_n = np.linalg.norm(rb, axis=1)
        cross = np.cross(ra, rb)
        denom = ra_n * rb_n * (ra_n * rb_n + np.einsum("ij,ij->i", ra, rb))
        # Avoid division by zero (collinear & outside the segment): clamp
        safe = denom > 1e-30
        scale = np.zeros_like(denom)
        scale[safe] = 1.0 / denom[safe]
        B += cross * scale[:, None] * (wire.current / (4 * np.pi))
    return B


def A_from_wire(rs: np.ndarray, wire: WirePath) -> np.ndarray:
    """Cartesian vector potential of a polyline filament.

    For a straight segment from a to b carrying current I from a to b,
    parameterise r' = a + u t_hat, u in [0, L], t_hat = (b-a)/L.
    Then
        A(r) = (mu_0 I / 4 pi) int_0^L du t_hat / |r - r'(u)|
             = (mu_0 I / 4 pi) t_hat * log( (|r-b| - s_b) / (|r-a| - s_a) )
    where s_a = t_hat . (r - a), s_b = t_hat . (r - b).
    (Sign convention checked: int du/sqrt((u-m)^2 + rho^2) from 0 to L
     = log[(L-m+|r-b|)/(-m+|r-a|)] = log[(|r-b|-s_b)/(|r-a|-s_a)] since
     -s_b = L - m and -s_a = -m.)
    """
    rs = np.atleast_2d(rs)
    A = np.zeros_like(rs)
    starts, ends, _ = wire_segments(wire)
    for a, b in zip(starts, ends):
        t = b - a
        L = np.linalg.norm(t)
        if L < 1e-15:
            continue
        t_hat = t / L
        ra = rs - a[None, :]
        rb = rs - b[None, :]
        ra_n = np.linalg.norm(ra, axis=1)
        rb_n = np.linalg.norm(rb, axis=1)
        sa = np.einsum("ij,j->i", ra, t_hat)
        sb = np.einsum("ij,j->i", rb, t_hat)
        num = rb_n - sb
        den = ra_n - sa
        # Near the wire itself num and den can be tiny; regularise using
        # the wire radius so A stays finite at the wire surface.
        floor = wire.radius
        num = np.maximum(num, floor)
        den = np.maximum(den, floor)
        logterm = np.log(num / den) * (wire.current / (4 * np.pi))
        A += np.outer(logterm, t_hat)
    return A


def wire_self_inductance(wire: WirePath) -> float:
    """Sum the standard self-inductance formula
        L_seg = (mu_0/2 pi) [ L (ln(2L/r) - 1) ]
    over straight segments, plus mutual inductance between non-adjacent
    segments via the Neumann formula (centroid quadrature).  This is
    approximate but adequate when segment lengths >> r_wire and
    segment-to-segment distances are comparable to segment lengths.

    mu_0 = 1 in our units.
    """
    starts, ends, dl = wire_segments(wire)
    n = starts.shape[0]
    L_arr = np.linalg.norm(dl, axis=1)
    r = wire.radius
    L_self = 0.0
    for i in range(n):
        Li = L_arr[i]
        if Li < 1e-12:
            continue
        L_self += (1.0 / (2 * np.pi)) * Li * (np.log(2 * Li / r) - 1.0)
    # Mutual inductance between non-adjacent pairs of segments via Neumann
    # formula approximated by centroid quadrature:
    #   M_ij = (mu_0/4pi) (dl_i . dl_j) / |c_i - c_j|
    centroids = 0.5 * (starts + ends)
    for i in range(n):
        for j in range(i + 2, n):  # skip self and immediate neighbour
            d = np.linalg.norm(centroids[i] - centroids[j])
            if d < 1e-12:
                continue
            M = (1.0 / (4 * np.pi)) * np.dot(dl[i], dl[j]) / d
            L_self += 2.0 * M
    return L_self


# ---------------------------------------------------------------------------
# Mutual-inductance kernel
# ---------------------------------------------------------------------------

def self_panel_factor(area: float, dchi: float, dt: float, s: float) -> float:
    """Approximate
        I_self = int_E int_E dS dS' / |r - r'|
    for a small rectangular panel of physical sides (s * dchi) x dt.
    Treats the panel as flat with sides a = s*dchi, b = dt.

    Analytic formula (Gradshteyn-Ryzhik 4.638):
      I = (1/3) * [ a^2 b * ln((b + sqrt(a^2+b^2))/a)
                  + a b^2 * ln((a + sqrt(a^2+b^2))/b)
                  + (a^3 + b^3 - (a^2+b^2)^{3/2}) / ... ]
    For our purposes use the compact form
      I_self(a,b) = a b * [ asinh(a/b) + (b/a) asinh(b/a) ] - (1/3)(a^3 + b^3 - (a^2+b^2)^{3/2})
    Actually the cleanest exact result:
      int_0^a int_0^a int_0^b int_0^b dx dy dx' dy' / sqrt((x-x')^2+(y-y')^2)
    = ... lengthy.  We use a numerical lookup based on Hammer-Smith formula:
        I_self(a,b)/(a*b)^{3/2} = g(b/a)
    with g(1) = 4/3 (2*sqrt(2) - ln(1+sqrt(2)) - ... )  ~ 3.525 for a=b.

    Implementation: use the closed form
      I = (2 a^2 b/3) asinh(b/a) + (2 a b^2/3) asinh(a/b)
        + (2/3)(sqrt(a^2 + b^2) (a^2 + b^2) - a^3 - b^3) / ... err.

    Rather than risk an algebra slip we numerically integrate once with
    Gauss-Legendre when needed, then cache.  However, for piecewise
    constant K with element area A, the self contribution to the
    integral int K.K /|r-r'| is just |K|^2 * I_self(a,b).
    """
    a = s * dchi
    b = dt
    # Use the closed-form expression for the four-fold integral
    # int_0^a int_0^a int_0^b int_0^b dx dx' dy dy' / sqrt((x-x')^2+(y-y')^2)
    # = (2/3) [ a^2 b sinh^{-1}(b/a) + a b^2 sinh^{-1}(a/b) ]
    #   + (2/3) [ (a^2 + b^2)^{3/2} - a^3 - b^3 ]
    # (Sourced from standard references on rectangular panel self-inductance,
    # see e.g. Stratton "Electromagnetic Theory" or Knoepfel "Magnetic Fields".)
    if a <= 0 or b <= 0:
        return 0.0
    asnh_ba = np.arcsinh(b / a)
    asnh_ab = np.arcsinh(a / b)
    term1 = (2.0 / 3.0) * (a * a * b * asnh_ba + a * b * b * asnh_ab)
    r2 = a * a + b * b
    term2 = (2.0 / 3.0) * (r2 ** 1.5 - a ** 3 - b ** 3)
    return term1 + term2


def assemble_mutual_matrix(mesh: Mesh,
                            near_threshold: float = 3.0) -> np.ndarray:
    """Build the n_el x n_el matrix
        G[i, j] = (1/(4 pi)) * area_i * area_j / |c_i - c_j|     (i != j)
                = (1/(4 pi)) * self_panel_factor(...)             (i == j)
    Then mutual inductance of basis K_i with basis K_j is
        M_ij = (mu_0/(4pi)) sum over... wait.

    Sorry, more carefully: the magnetic energy of surface currents is
        W = (1/(8pi)) int int K(r).K(r')/|r-r'| dS dS'   (mu_0=1)
    With piecewise-constant K per element K_e, area A_e:
        W = (1/(8pi)) sum_{e,e'} (K_e . K_e') G(e, e')
    where
        G(e, e') = int_{E_e} int_{E_{e'}} dS dS' / |r-r'|.
    For e != e' and well-separated panels, G(e,e') ~ A_e A_e' / |c_e - c_e'|
    (1-point centroid quadrature).  For close panels use higher-order.
    For self e == e' use self_panel_factor.

    We return G as a dense (n_el, n_el) array.  Memory: 8 * n_el^2 bytes.
    For n_el = 1500 that's 18 MB; for n_el = 6400 that's 330 MB.
    """
    n = mesh.elem_xyz.shape[0]
    G = np.empty((n, n), dtype=np.float64)
    centroids = mesh.elem_xyz
    areas = mesh.elem_area
    # Vectorised distances
    diff = centroids[:, None, :] - centroids[None, :, :]
    dist = np.sqrt(np.einsum("ijk,ijk->ij", diff, diff))
    # Off-diagonal: A_i A_j / dist
    np.fill_diagonal(dist, 1.0)  # avoid div by zero
    G = (areas[:, None] * areas[None, :]) / dist
    # Diagonal: self panel
    diag = np.array([self_panel_factor(mesh.elem_area[i],
                                       mesh.elem_dchi[i],
                                       mesh.elem_dt[i],
                                       mesh.elem_s[i]) for i in range(n)])
    np.fill_diagonal(G, diag)
    return G


def assemble_wire_coupling(mesh: Mesh, wire: WirePath) -> np.ndarray:
    """Vector b_e = int_{E_e} dS A_wire(r) . (something).

    The cross term in the energy is
        W_cross = int K(r) . A_wire(r) dS
    For piecewise constant K per element with K_e:
        W_cross = sum_e K_e . [ int_{E_e} A_wire(r) dS ]
                ~ sum_e K_e . A_wire(c_e) * area_e        (centroid quadrature)
    We return the (n_el, 3) array v[e] = A_wire(c_e) * area_e.
    """
    A_w = A_from_wire(mesh.elem_xyz, wire)
    return A_w * mesh.elem_area[:, None]


# ---------------------------------------------------------------------------
# Linear-system assembly via energy minimisation
# ---------------------------------------------------------------------------

def assemble_psi_linear_system(mesh: Mesh, wire: WirePath):
    """Set up M_nodal . psi_int = rhs.

    K (Cartesian, per element) is a linear function of nodal psi:
        K_e = sum_n  T[e, n, :] * psi_n
    where T[e, n, :] is a 3-vector contribution per node-per-element.

    The non-trivial T entries are only the 4 corner nodes of element e:
        T[e, corners[e, k], :] = grad_chi_coef[e, k] * K_chi_coef_cart[e]
                                + grad_t_coef[e, k]  * K_t_coef_cart[e]

    Sheath energy (mu_0=1):
        W_sheath = (1/(8 pi)) sum_{e, e'} (K_e . K_{e'}) G(e, e')
                 = (1/(8 pi)) psi^T A psi
    where A_nm = sum_{e, e'} G(e,e') * (T[e,n,:] . T[e',m,:]).

    Cross term with wire:
        W_cross = sum_e K_e . v[e]  (where v[e] = A_wire(c_e)*area_e)
                = sum_n psi_n * c_n,
        c_n = sum_e (T[e,n,:] . v[e]).

    Total quadratic problem in psi:
        W = (1/(8 pi)) psi^T A psi + psi^T c + L_wire_self
    delta W / delta psi_int = 0  =>  (1/(4 pi)) A_{int,int} psi_int
                                   = -A_{int,bnd} psi_bnd / (4 pi) - c_int

    Returns:
      M_int_int (n_int, n_int),
      rhs_int   (n_int,),
      psi_bnd_full (n_nodes,) with Dirichlet values inserted,
      A (full nodal matrix, before Dirichlet split),
      c (n_nodes,)
    """
    n_el = mesh.elem_xyz.shape[0]
    n_n = mesh.n_nodes

    # Build T as sparse representation: for each element, 4 (node, K-vec) pairs.
    # T_eN[e, k, :] = corner k's contribution to K_e (Cartesian).
    T_eK = (mesh.grad_chi_coef[:, :, None] * mesh.K_chi_coef_cart[:, None, :]
            + mesh.grad_t_coef[:, :, None] * mesh.K_t_coef_cart[:, None, :])
    # shape: (n_el, 4, 3)
    corners = mesh.elem_node_idx           # (n_el, 4)

    G = assemble_mutual_matrix(mesh)       # (n_el, n_el)
    v = assemble_wire_coupling(mesh, wire) # (n_el, 3)

    # A_nm  via efficient accumulation.  We'll do it as
    #   A = sum_{ne, e'} sum_{k,l} T_eK[e,k,:].T_eK[e',l,:] G[e,e']
    #       indexed by (corners[e,k], corners[e',l]).
    # This is O(n_el^2 * 16) operations.  For n_el=1500: 36M ops, ok.
    A_full = np.zeros((n_n, n_n), dtype=np.float64)
    c_full = np.zeros(n_n, dtype=np.float64)

    # Cross-term c_n = sum_e sum_k delta_{corners[e,k]=n} T_eK[e,k,:].v[e]
    Tdotv = np.einsum("ekc,ec->ek", T_eK, v)     # (n_el, 4)
    np.add.at(c_full, corners.flatten(), Tdotv.flatten())

    # A matrix: sum_{e,e'} K_e . K_e' G[e,e']
    # Decompose: K_e = sum_k T[e,k] psi[corners[e,k]]
    # Define D[e, n, :] = sum_k delta_{corners[e,k]=n} T[e,k,:]
    # But that's still n_el x n_n x 3.  Better: directly accumulate.
    # Block: for each e, build B[e, :, :] = T_eK[e] @ ones with corner labels.
    # Then A += B[e].T @ (G[e, e'] * B[e']) but that's basically O(n_el^2 n_n).
    #
    # Cleaner: form K_basis matrix B of shape (n_el, 3, n_n) sparsely.
    # We'll do it by computing the full (n_el, 3) "synthesised K" for each
    # psi unit vector.  That's n_n forward solves: too expensive.
    #
    # Vectorised path: form a (n_el, n_n, 3) sparse tensor S where
    # S[e, n, :] = sum_k delta_{corners[e,k]=n} T_eK[e, k, :].
    # Then A = sum_{e,e'} G[e,e'] * (S[e,:,:] @ S[e',:,:].T)_summed_over_coords
    # Practically, store S as scipy.sparse in shape (3*n_el, n_n)?  Yes:
    # Stack S over axis 1 to get S_flat[e*3 + j, n], then A = S_flat.T @ G_kron3 @ S_flat?
    # But G is block-shared.  Use a direct loop with vectorisation:
    #
    # For each n', do:
    #   K_n'(e) = T[e,:,:].psi_n'_at_corners(e)    # (3,) per e
    #   y(e)    = sum_{e'} G[e,e'] * K_n'(e')      # (3,) per e
    #   A[:, n'] is then n -> sum_e T[e,:,n] . y[e]
    # Cost: O(n_n * n_el^2).  For n_n ~ 6000, n_el ~ 1500, that's 5e10 ops.
    # Too slow.
    #
    # Better: collapse element K-vectors to a vector of size 3*n_el, define
    # the kernel matrix G3 of shape (3*n_el, 3*n_el) = block-diagonal with
    # G[e,e'] * I_3 in each (e,e') block, then
    #   K_3 = S @ psi,  S in R^{3 n_el x n_n} sparse,
    #   A = S.T @ G3 @ S, all linear operations.
    # G3 is (3 n_el)^2 dense memory: for n_el=1500, 4500^2 = 20M entries
    # = 160 MB.  Acceptable.
    # For n_el=4000, that's 1.5 GB.  Watch.
    # A_full[n, m] = sum_j sum_{e, e'} S_j[e, n] G[e, e'] S_j[e', m]
    # where S_j is the sparse n_el x n_n matrix with values T_eK[e, k, j] at
    # (e, corners[e, k]).  We exploit block-diagonal structure of G3 to avoid
    # ever materialising the 3*n_el x 3*n_el matrix.
    from scipy.sparse import csr_matrix
    A_full = np.zeros((n_n, n_n), dtype=np.float64)
    e_idx = np.repeat(np.arange(n_el), 4)
    col_idx = corners.flatten()
    for j in range(3):
        vals = T_eK[:, :, j].flatten()
        S_j = csr_matrix((vals, (e_idx, col_idx)), shape=(n_el, n_n))
        # G @ S_j is dense (n_el, n_n); S_j.T @ (G @ S_j) is dense (n_n, n_n)
        GSj = G @ S_j           # csr_matrix dense result is np.ndarray
        if not isinstance(GSj, np.ndarray):
            GSj = np.asarray(GSj.todense())
        A_full += np.asarray(S_j.T @ GSj)

    # Apply 1/(4 pi) prefactor for delta W = 0 condition (after symm in 1/(8pi)*2).
    # delta W / delta psi_n = (1/(4pi)) sum_m A[n,m] psi_m + c[n] = 0
    A_full = A_full / (4 * np.pi)
    # c_full picks up no extra 1/(4 pi) since we already absorbed it in
    # A_wire (which itself has the mu_0/(4 pi) factor).  Verify:
    # A_wire(r) = (mu_0/(4 pi)) int dl'/|r-r'|, mu_0=1.  Yes.
    # So c_full is already in the units needed.

    # Dirichlet split.
    psi_bnd_full = np.zeros(n_n)
    psi_bnd_full[mesh.right_node_idx] = +0.5
    psi_bnd_full[mesh.left_node_idx] = -0.5

    int_idx = mesh.interior_node_idx
    bnd_idx = np.setdiff1d(np.arange(n_n), int_idx)

    A_ii = A_full[np.ix_(int_idx, int_idx)]
    A_ib = A_full[np.ix_(int_idx, bnd_idx)]
    rhs = -A_ib @ psi_bnd_full[bnd_idx] - c_full[int_idx]
    return A_ii, rhs, psi_bnd_full, A_full, c_full


# ---------------------------------------------------------------------------
# Solve and post-process
# ---------------------------------------------------------------------------

def solve_psi(mesh: Mesh, wire: WirePath, verbose: bool = False):
    t0 = time.time()
    A_ii, rhs, psi_bnd_full, A_full, c_full = assemble_psi_linear_system(mesh, wire)
    t1 = time.time()
    if verbose:
        print(f"Assembly: {t1 - t0:.2f} s, A_ii shape {A_ii.shape}", flush=True)
    # Solve
    psi_int = la_solve(A_ii, rhs, assume_a="sym")
    t2 = time.time()
    if verbose:
        print(f"Solve:    {t2 - t1:.2f} s", flush=True)
    psi = psi_bnd_full.copy()
    psi[mesh.interior_node_idx] = psi_int
    return psi, dict(A_full=A_full, c_full=c_full, assemble_time=t1 - t0,
                     solve_time=t2 - t1)


def compute_Lp(psi: np.ndarray, mesh: Mesh, wire: WirePath,
               full_info: dict) -> float:
    """L_p = 2 * W_total = (1/(4 pi)) psi^T A psi + 2 psi^T c + L_wire_self.

    Derivation: total magnetic energy
       W = (1/(8 pi)) psi^T A_quad psi + psi^T c + (1/2) L_wire_self
    where A_quad = (4 pi) A_full as we stored A_full = A_quad/(4 pi).
       W = (1/2) psi^T A_full psi + psi^T c + (1/2) L_wire_self
    L_p (mu_0 = 1, unit current) is L = 2 W:
       L_p = psi^T A_full psi + 2 psi^T c + L_wire_self.
    """
    A_full = full_info["A_full"]
    c_full = full_info["c_full"]
    L_sheath = psi @ (A_full @ psi)
    L_cross = 2.0 * (psi @ c_full)
    L_wire = wire_self_inductance(wire)
    return L_sheath + L_cross + L_wire


def compute_Veff(psi: np.ndarray, mesh: Mesh, wire: WirePath,
                 n_s: int = 20, n_phi: int = 50, n_z: int = 20) -> float:
    """Volume overlap
        V_eff = (R/mu_0) int_{V_mag} A_phi(s, phi, z) * s ds dphi dz
    with R = mu_0 = 1.  Integrand:
        (1/s) * A_phi   actually no.  B_0 = (R/s) phi_hat, so
        V_eff = (1/(mu_0 B_max)) int B_0 . A1 dV
              = (1/mu_0) int (R/s) A1_phi * s ds dphi dz
              = (R/mu_0) int A1_phi ds dphi dz
        no s factor.

    With R = 1, mu_0 = 1:
        V_eff = int A1_phi(s, phi, z) ds dphi dz   over V_mag.

    A1 = A_sheath + A_wire, where
        A_sheath(r) = (1/(4 pi)) sum_e K_e * area_e / |r - c_e|
    and A_wire by finite-segment formulas.
    """
    alpha, beta = mesh.alpha, mesh.beta
    # Gauss-Legendre nodes in [0, 1]
    s_nodes, s_w = np.polynomial.legendre.leggauss(n_s)
    z_nodes, z_w = np.polynomial.legendre.leggauss(n_z)
    phi_nodes, phi_w = np.polynomial.legendre.leggauss(n_phi)
    # Map to physical ranges
    s_vals = 1.0 + 0.5 * beta * (s_nodes + 1)
    s_jac = 0.5 * beta
    z_vals = 0.5 * alpha * (z_nodes + 1)
    z_jac = 0.5 * alpha
    phi_vals = np.pi * (phi_nodes + 1)
    phi_jac = np.pi

    # Build grid of points and integrate.
    K = K_from_psi(psi, mesh)         # (n_el, 3)
    centroids = mesh.elem_xyz         # (n_el, 3)
    areas = mesh.elem_area            # (n_el,)
    weights_K = areas / (4 * np.pi)   # (n_el,) -- (mu_0=1)/4 pi * area

    total = 0.0
    # Loop over phi to keep memory under control.  At each phi compute the
    # 2D (s, z) integral.
    for ip, (phi_v, w_phi) in enumerate(zip(phi_vals, phi_w)):
        # 2D (s, z) grid at this phi
        ss, zz = np.meshgrid(s_vals, z_vals, indexing="ij")
        ws, wz = np.meshgrid(s_w, z_w, indexing="ij")
        xs = ss * np.cos(phi_v)
        ys = ss * np.sin(phi_v)
        zs = zz
        pts = np.stack([xs.flatten(), ys.flatten(), zs.flatten()], axis=1)
        # n_pts x 3
        # A_sheath(p) . phi_hat
        # A_sheath = sum_e K_e * area_e/(4 pi |p - c_e|)
        diff = pts[:, None, :] - centroids[None, :, :]
        dist = np.sqrt(np.einsum("ijk,ijk->ij", diff, diff))
        # A_x and A_y
        Ax = (K[None, :, 0] * weights_K[None, :] / dist).sum(axis=1)
        Ay = (K[None, :, 1] * weights_K[None, :] / dist).sum(axis=1)
        # A_wire at these points
        Aw = A_from_wire(pts, wire)
        Ax += Aw[:, 0]
        Ay += Aw[:, 1]
        # phi_hat at azimuth phi_v: (-sin phi, cos phi, 0).
        A_phi = -Ax * np.sin(phi_v) + Ay * np.cos(phi_v)
        # Multiply by quadrature weights and Jacobian.
        # int_{V_mag} A_phi ds dphi dz = sum w_s w_z w_phi A_phi(s,phi,z) * s_jac z_jac phi_jac
        local = (A_phi.reshape(s_vals.size, z_vals.size)
                 * ws * wz).sum()
        total += w_phi * local
    total *= s_jac * z_jac * phi_jac
    return total


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

@dataclass
class RunResult:
    alpha: float
    beta: float
    w_nominal: float
    w_eff: float
    n_elements: int
    F: float
    ell: float
    L_p: float
    V_eff: float
    matrix_cond: float
    wall_time_s: float


def estimate_w_eff(mesh: Mesh) -> float:
    """Effective slit width = element size of the elements adjacent to chi=0
    or chi=Phi.  Return the larger of the two for conservativeness."""
    h_lo = mesh.chi_nodes[1] - mesh.chi_nodes[0]
    h_hi = mesh.chi_nodes[-1] - mesh.chi_nodes[-2]
    return max(h_lo, h_hi)


def run_one(alpha: float,
            beta: float,
            w: float = 0.03,
            n_t_per_seg: int = 20,
            chi_nodes_override: Optional[np.ndarray] = None,
            wire_type: str = "production",
            verbose: bool = False) -> Tuple[RunResult, dict]:
    t0 = time.time()
    mesh = build_mesh(alpha, beta, w=w, n_t_per_seg=n_t_per_seg,
                      chi_nodes_override=chi_nodes_override)
    if wire_type == "production":
        wire = production_wire(alpha, w=w)
    elif wire_type == "linked":
        wire = linked_wire(alpha, w=w)
    elif wire_type == "none":
        wire = empty_wire()
    else:
        raise ValueError(wire_type)
    psi, info = solve_psi(mesh, wire, verbose=verbose)
    L_p = compute_Lp(psi, mesh, wire, info)
    V_eff = compute_Veff(psi, mesh, wire)
    t1 = time.time()
    cond = np.linalg.cond(info["A_full"][np.ix_(mesh.interior_node_idx,
                                                mesh.interior_node_idx)])
    res = RunResult(
        alpha=alpha, beta=beta,
        w_nominal=w, w_eff=estimate_w_eff(mesh),
        n_elements=mesh.elem_xyz.shape[0],
        F=V_eff, ell=L_p,
        L_p=L_p, V_eff=V_eff,
        matrix_cond=cond,
        wall_time_s=t1 - t0,
    )
    if verbose:
        print(f"alpha={alpha}, beta={beta}: F={V_eff:.4e}, ell={L_p:.4e},"
              f" n_el={res.n_elements}, t={res.wall_time_s:.1f}s,"
              f" cond={cond:.2e}", flush=True)
    return res, {"psi": psi, "mesh": mesh, "info": info}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha", type=float, required=True)
    ap.add_argument("--beta", type=float, required=True)
    ap.add_argument("--w", type=float, default=0.03)
    ap.add_argument("--nt", type=int, default=20, help="n_t_per_seg")
    ap.add_argument("--out", type=str, required=True)
    ap.add_argument("--wire", type=str, default="production",
                     choices=["production", "linked", "none"])
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    res, _ = run_one(args.alpha, args.beta, w=args.w,
                     n_t_per_seg=args.nt, wire_type=args.wire,
                     verbose=args.verbose)
    out = asdict(res)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()

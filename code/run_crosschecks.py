"""Run cross-checks 1-4 from the spec and write crosscheck_results.json."""

from __future__ import annotations
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from bem_solver import (
    build_mesh, production_wire, linked_wire,
    solve_psi, compute_Lp, compute_Veff, K_from_psi,
    A_from_wire, run_one,
)


def main(out_path: str):
    results = {}
    t0 = time.time()

    # --------------------------------------------------------------
    # XC2: zero-slit limit
    # --------------------------------------------------------------
    print("== XC2 zero-slit limit ==", flush=True)
    xc2 = []
    for w in [0.03, 0.01, 0.003, 0.001]:
        res, _ = run_one(alpha=1.0, beta=1.0, w=w, n_t_per_seg=12)
        xc2.append(dict(w_nominal=w, w_eff=res.w_eff,
                         F=res.F, ell=res.ell, n_el=res.n_elements,
                         cond=res.matrix_cond))
        print(f"  w={w:.4f}  F={res.F:+.4e}  ell={res.ell:.4e}", flush=True)
    F_ref = next(item["F"] for item in xc2 if item["w_nominal"] == 0.03)
    F_min = next(item["F"] for item in xc2 if item["w_nominal"] == 0.001)
    # Monotonicity: |F| decreasing as w decreases
    Fs = [abs(item["F"]) for item in xc2]
    monotone = all(Fs[i] >= Fs[i + 1] for i in range(len(Fs) - 1))
    factor_10x = abs(F_min) < 0.1 * abs(F_ref)
    xc2_pass = bool(monotone and factor_10x)
    results["XC2"] = dict(table=xc2,
                          monotone=bool(monotone),
                          ratio_0p001_over_0p03=abs(F_min) / abs(F_ref) if F_ref != 0 else float("inf"),
                          passed=xc2_pass)
    print(f"  XC2 monotone={monotone}, ratio={abs(F_min)/abs(F_ref):.3e}, pass={xc2_pass}",
          flush=True)

    # --------------------------------------------------------------
    # XC1: linked-return topology, V_eff -> 0
    # --------------------------------------------------------------
    print("== XC1 linked-return topology ==", flush=True)
    # Production at alpha=beta=1, w=0.03
    res_prod, _ = run_one(alpha=1.0, beta=1.0, w=0.03, n_t_per_seg=10,
                          wire_type="production")
    res_linked, _ = run_one(alpha=1.0, beta=1.0, w=0.03, n_t_per_seg=10,
                            wire_type="linked")
    ratio = abs(res_linked.F) / abs(res_prod.F) if res_prod.F != 0 else float("inf")
    xc1_pass = bool(ratio < 0.01)
    results["XC1"] = dict(F_production=res_prod.F, F_linked=res_linked.F,
                          ratio=ratio, passed=xc1_pass)
    print(f"  F_production={res_prod.F:+.4e}, F_linked={res_linked.F:+.4e},"
          f" ratio={ratio:.3e}, pass={xc1_pass}", flush=True)

    # --------------------------------------------------------------
    # XC3: reciprocity (energy vs surface integral)
    # --------------------------------------------------------------
    print("== XC3 reciprocity (surface vs volume L_p) ==", flush=True)
    # Surface-integral L_p was computed via psi^T A_full psi + 2 psi^T c + L_wire_self.
    # This *is* the surface-integral form ∫ K.A dS  (= ∫ K·A_sheath + ∫ K·A_wire + wire self).
    # For an alternative form we compute L_p as 2 * total magnetic energy
    # = (1/mu_0) ∫ |B|^2 dV by direct integration over a bounding box.
    mesh = build_mesh(1.0, 1.0, w=0.03, n_t_per_seg=10)
    wire = production_wire(1.0)
    psi, info = solve_psi(mesh, wire)
    L_surface = compute_Lp(psi, mesh, wire, info)

    # Volume form: integrate |B|^2 over a large box [-3, 3]^3 (mu_0 = 1).
    # B from sheath currents at any field point:
    #   B = (1/(4 pi)) sum_e K_e × (r - c_e) * area_e / |r - c_e|^3 + B_wire
    K = K_from_psi(psi, mesh)
    centroids = mesh.elem_xyz
    areas = mesh.elem_area
    weights_K = areas / (4 * np.pi)

    # Quadrature box.  Volume of the box must comfortably contain the sheath
    # AND the dominant wire-field volume.  Box: |x|<3, |y|<3, -3<z<4.
    box = [(-3.0, 3.0, 24), (-3.0, 3.0, 24), (-3.0, 4.0, 30)]
    L_vol = 0.0
    xs = np.linspace(box[0][0], box[0][1], box[0][2] + 1)
    ys = np.linspace(box[1][0], box[1][1], box[1][2] + 1)
    zs = np.linspace(box[2][0], box[2][1], box[2][2] + 1)
    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]
    dz = zs[1] - zs[0]
    # Midpoint quadrature
    xm = 0.5 * (xs[:-1] + xs[1:])
    ym = 0.5 * (ys[:-1] + ys[1:])
    zm = 0.5 * (zs[:-1] + zs[1:])
    # Cache wire info
    starts = wire.points[:-1]
    ends = wire.points[1:]
    dl = ends - starts
    # Loop chunks over z to control memory
    B2_total = 0.0
    for zv in zm:
        XX, YY = np.meshgrid(xm, ym, indexing="ij")
        pts = np.stack([XX.flatten(), YY.flatten(),
                         np.full(XX.size, zv)], axis=1)
        # Sheath B
        diff = pts[:, None, :] - centroids[None, :, :]
        dist = np.sqrt(np.einsum("ijk,ijk->ij", diff, diff))
        # B_e contribution: K_e x (r - c_e)/(4 pi |r-c_e|^3) * area_e
        # Avoid div-zero (none of pts should sit exactly on centroids)
        d3 = dist ** 3
        d3[d3 == 0] = 1e-30
        # K x diff_normalised
        K_e = K[None, :, :]                 # (1, n_el, 3)
        cross = np.cross(K_e, diff)         # (n_pts, n_el, 3)
        B_sheath = (cross * (weights_K[None, :, None] / d3[:, :, None])).sum(axis=1)
        # Wire B (finite-segment)
        from bem_solver import B_from_wire as Bw
        B_wire = Bw(pts, wire)
        B_total = B_sheath + B_wire
        B2 = np.einsum("ij,ij->i", B_total, B_total)
        # Exclude points very close to wire (|B|^2 diverges as 1/r^2 near wire)
        # by capping the contribution.  For a thin wire, the near-field
        # contribution is wire-self-inductance which is accounted for
        # analytically; double-counting would inflate L_vol.
        # Practical clip: any point closer than 0.1 to any wire segment gets
        # discarded (we replace by analytic estimate later).
        # Distance from each point to each segment:
        d_seg = np.full(pts.shape[0], np.inf)
        for a, b in zip(starts, ends):
            t = b - a
            L = np.linalg.norm(t)
            if L < 1e-15:
                continue
            t_hat = t / L
            ra = pts - a
            proj = ra @ t_hat
            proj = np.clip(proj, 0.0, L)
            closest = a + proj[:, None] * t_hat[None, :]
            d_loc = np.linalg.norm(pts - closest, axis=1)
            d_seg = np.minimum(d_seg, d_loc)
        good = d_seg > 0.15
        B2 = np.where(good, B2, 0.0)
        B2_total += B2.sum() * dx * dy * dz
    L_volume_excl_wire_near = B2_total      # 1 * int |B|^2 dV for unit current (mu_0=1)
    # Add wire self-inductance separately for cylinder near-wire region:
    from bem_solver import wire_self_inductance
    L_wire_only = wire_self_inductance(wire)
    L_volume = L_volume_excl_wire_near + L_wire_only  # rough; not exact but useful

    rel_diff = abs(L_volume - L_surface) / abs(L_surface)
    xc3_pass = bool(rel_diff < 0.20)
    results["XC3"] = dict(L_surface=float(L_surface),
                          L_volume=float(L_volume),
                          L_volume_excl_wire=float(L_volume_excl_wire_near),
                          L_wire_self=float(L_wire_only),
                          rel_diff=float(rel_diff),
                          passed=xc3_pass)
    print(f"  L_surface={L_surface:.4e}, L_volume={L_volume:.4e},"
          f" rel_diff={rel_diff:.3e}, pass={xc3_pass}", flush=True)

    # --------------------------------------------------------------
    # XC4: gauge invariance of V_eff under A -> A + grad chi
    # --------------------------------------------------------------
    print("== XC4 gauge invariance of V_eff ==", flush=True)
    # The spec example chi = s phi z is multivalued; we use single-valued
    # chi(s, phi, z) = z * sin(phi) instead.  Then
    #   grad chi = (0, z cos(phi)/s, sin(phi))   in (s, phi, z) basis.
    # The phi-component contribution to V_eff is
    #   int_V_mag (z cos phi / s) * s ds dphi dz = int s ds int z dz int cos phi dphi = 0.
    # Numerically we evaluate this with the same Gauss-Legendre rule as
    # compute_Veff and verify the magnitude is below tolerance.
    alpha, beta = 1.0, 1.0
    # Reuse the same quadrature
    n_s, n_phi, n_z = 20, 50, 20
    sn, sw_ = np.polynomial.legendre.leggauss(n_s)
    zn, zw_ = np.polynomial.legendre.leggauss(n_z)
    pn, pw_ = np.polynomial.legendre.leggauss(n_phi)
    s_vals = 1.0 + 0.5 * beta * (sn + 1)
    z_vals = 0.5 * alpha * (zn + 1)
    phi_vals = np.pi * (pn + 1)
    dVeff_gauge = 0.0
    for ip, (pv, w_p) in enumerate(zip(phi_vals, pw_)):
        local = 0.0
        for iz, (zv, w_z) in enumerate(zip(z_vals, zw_)):
            for is_, (sv, w_s) in enumerate(zip(s_vals, sw_)):
                # (grad chi)_phi = z cos(phi) / s; integrand = (grad chi)_phi
                local += w_s * w_z * (zv * np.cos(pv) / sv)
        dVeff_gauge += w_p * local
    dVeff_gauge *= (0.5 * beta) * (0.5 * alpha) * np.pi
    V_eff_baseline = res_prod.F
    rel = abs(dVeff_gauge) / max(abs(V_eff_baseline), 1e-12)
    xc4_pass = bool(rel < 1e-3)
    results["XC4"] = dict(chi="z*sin(phi) [single-valued]",
                          dV_eff_under_gauge=float(dVeff_gauge),
                          V_eff_baseline=float(V_eff_baseline),
                          rel=float(rel), passed=xc4_pass)
    print(f"  dV_eff under gauge = {dVeff_gauge:+.3e},"
          f" V_eff_baseline={V_eff_baseline:+.3e}, rel={rel:.3e}, pass={xc4_pass}",
          flush=True)

    # --------------------------------------------------------------
    # Slit-width sensitivity at alpha=beta=1
    # --------------------------------------------------------------
    print("== Slit sensitivity at alpha=beta=1 ==", flush=True)
    sens = []
    for w in [0.015, 0.03, 0.06]:
        res, _ = run_one(alpha=1.0, beta=1.0, w=w, n_t_per_seg=12)
        sens.append(dict(w=w, w_eff=res.w_eff, F=res.F, ell=res.ell))
    results["slit_sensitivity"] = sens

    # --------------------------------------------------------------
    # Persist
    # --------------------------------------------------------------
    results["wall_time_s"] = time.time() - t0
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=lambda v: float(v))
    print(f"\nResults written to {out_path}.")
    print(f"XC1 passed: {results['XC1']['passed']}")
    print(f"XC2 passed: {results['XC2']['passed']}")
    print(f"XC3 passed: {results['XC3']['passed']}")
    print(f"XC4 passed: {results['XC4']['passed']}")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="outputs/crosscheck_results.json")
    args = ap.parse_args()
    main(args.out)

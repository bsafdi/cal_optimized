# Patch: remove external return wire from BEM calculation

The v2 spec's external wire path created a parasitic inductance floor (ℓ_0 ≈ 19
μ_0R) that dominated L_p and forced V_eff to be a residual on a ~99.9% cancellation,
which sat below the kernel quadrature noise floor. We're now isolating the sheath
in vacuum so the decision criterion can bite on the sheath physics alone.

This is a small patch to the existing `bem_solver.py` (current commit: e70bac2,
the post-sign-fix version). Mesh generation, basis functions, BEM matrix assembly
on the sheath, slit-edge node identification, the four-segment slit topology, and
the Δψ = 1 Dirichlet constraint are all unchanged.

## Changes to `bem_solver.py`

1. **Remove wire from the linear-system RHS.** The BEM matrix currently has
   `RHS[i] = -n̂(r_i) · B_wire(r_i)` on each sheath element centroid. Set that
   contribution to zero. The Δψ = 1 Dirichlet rows are unchanged. The system
   becomes:
   - Meissner BCs: n̂ · B_K = 0 on each sheath element (homogeneous in K)
   - Slit constraint: Δψ = 1 across the slit (inhomogeneous)
   The solution is the unique minimum-energy K consistent with both constraints.

2. **Remove wire from L_p.** Currently L_p has two pieces:
   - Surface integral `∫_S K · A dS` where `A` includes wire contributions
   - Additive `L_wire,self` term
   Two edits:
   - In computing `A` on the sheath for the surface integral, use only the
     sheath's own surface currents (no wire contribution). Equivalent compact
     form:
     ```
     L_p = (μ_0/4π) ∫∫_S [K(r)·K(r')] / |r - r'| dS dS'
     ```
   - Delete the `L_wire,self` additive term entirely.

3. **Remove wire from V_eff.** Currently `V_eff = (R/μ_0) ∫_{V_mag} A_φ s ds dφ dz`
   with `A_φ` from both sheath and wire. Use `A_φ` from sheath surface currents
   only (Biot-Savart on K). No other change to the volume quadrature.

4. **No other code changes.** Keep the adaptive φ-mesh, the near-singular kernel
   treatment for adjacent panels, the slit Dirichlet handling, and the four-segment
   slit topology exactly as in e70bac2. The sign-convention fix from that commit
   stays.

## Cross-checks to re-run

- **XC1 (linked-return topology):** SKIP. With no wire, the linked-return
  topology is not realised. Mark as not applicable in `crosscheck_results.json`.
- **XC2 (zero-slit limit):** RE-RUN at α = β = 1, sweeping w ∈ {0.03, 0.01, 0.003,
  0.001}. Pass criterion unchanged: monotone decrease, |V_eff(w=0.001)| <
  0.1 × |V_eff(w=0.03)|.
- **XC3 (reciprocity):** RE-RUN. Surface integral form `(μ_0/4π) ∫∫ K·K' /|r-r'|`
  vs. volume integral form `(1/μ_0) ∫ |B|² dV` over a large box (now no wire
  self-energy to subtract). Pass criterion: < 5% relative difference.
- **XC4 (gauge invariance):** RE-RUN. Add `∇χ` with `χ = z sin φ` to the computed
  `A` and verify `V_eff` unchanged. Pass criterion: relative shift < 10⁻⁴.

All three (XC2, XC3, XC4) must pass before the production sweep.

## Production sweep

After cross-checks pass, re-run the full 18-point sweep:
- α ∈ {0.5, 1.0, 2.0}
- β ∈ {0.20, 0.12, 0.08, 0.05, 0.035, 0.025}
- w = 0.03 (nominal; record w_eff per run)

For each α, refit:
```
F(β) = f_1 β + f_2 β^2 + f_3 β^3
ℓ(β) = ℓ_0 + ℓ_1 β + ℓ_2 β^2
```
Bootstrap residuals (2000 samples) for uncertainties. Apply the v2 decision
criterion unchanged:
- f_1 statistically nonzero if |f_1| > 3 × σ_boot(f_1)
- ℓ_0 consistent with zero if RSS(constrained ℓ_0=0) ≤ 2 × RSS(unconstrained)

## Expected behavior (so you can flag deviations)

These are predictions, not requirements. Report what you find regardless.

- **L_p magnitudes:** ℓ should drop from ~19 to something of order unity or less.
  The wire floor is gone; what remains is the sheath's intrinsic self-inductance.
- **L_p β-scaling:** if the sheath alone has favorable thin-annulus scaling, ℓ_0
  should now be small (consistent with zero) and ℓ_1 β should dominate. If
  there's an intrinsic sheath inductance floor (e.g., from slit-edge or finite
  cap effects), ℓ_0 will still be nonzero even without the wire — that would be
  a real physics finding, not a spec error.
- **V_eff:** the ~99.9% cancellation between wire A_φ and induced sheath A_φ is
  gone, so V_eff should resolve well above the kernel noise floor. F(β) should
  now show clean β-dependence (monotone, fittable). If F still doesn't resolve,
  the kernel quadrature needs upgrading and that's a separate task.
- **Decision per α:** for the first time this is meaningful. Apply the criterion
  honestly. If results differ between α values, report the α-dependence.

## Outputs

Update in `/mnt/user-data/outputs/`:
- `production_sweep.csv` — new F, ℓ values for the 18 (α, β) points
- `fit_results.json` — new f_1, ℓ_0, ℓ_1 with bootstrap uncertainties, new
  decision per α
- `crossover_analysis.json` — recompute crossover R_t* only if at least one α
  passes the decision criterion. Otherwise mark "not applicable" with reason.
- `crosscheck_results.json` — XC2, XC3, XC4 results; XC1 marked as N/A.
- `report.md` — updated executive summary stating:
  - The wire was removed; this run is the sheath-only result
  - Pass/fail of XC2, XC3, XC4
  - Fit coefficients per α with uncertainties
  - Decision per α, with reasoning
  - Crossover R_t* (if applicable) in physical units for V_mag ∈ {0.1, 1, 10} m³

## Stop conditions

Don't tune. If you hit one of these, stop and report:

- Cross-check XC2, XC3, or XC4 fails after the patch. Means the wire removal
  introduced a bug; debug the patch.
- F(β) does not converge with mesh refinement at α = β = 1 (test at
  nt ∈ {10, 15, 20, 25}). Means the kernel quadrature is still inadequate even
  without the wire; this is a real result that needs reporting, and the
  upgrade path (4-point Gauss + analytic near-singular split) is then warranted
  as a follow-up task.
- ℓ(β) is non-monotone or doesn't fit a polynomial. Means the sheath physics
  itself has structure we didn't anticipate; report the raw data.

In any of these cases, `report.md` should document exactly what was seen and
what is and isn't concluded. An honest "wire removed; sheath physics doesn't
resolve cleanly at current resolution" is more useful than a forced fit.

## Estimated time

Patch + cross-checks: ~1 hour. Production sweep: ~5 minutes (same compute as
before; the matrix size doesn't change, just the RHS and a few additive terms).
Refits + report: ~30 minutes. Total: under 2 hours of agent time.

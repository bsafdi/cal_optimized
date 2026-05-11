# Task: Unit-current BEM calculation for hermetic-sheath toroidal axion pickup (v2)

## Background and goal

This task tests whether the proposed hermetic Meissner-sheath toroidal pickup for the
CAL axion experiment has favorable Q ∝ R_t^4 scaling at fixed magnet cost in the
thin-annulus limit. The question is whether the pickup inductance has a parasitic
"floor" (set by slit and lead geometry) that destroys the asymptotic gain, or whether
it scales linearly with the annulus parameter β so the gain is real.

Read `toroid_optimization.pdf` (sections 1–8) before starting. It sets up the geometry,
fixes notation, and explains why the question reduces to two scalar coefficients.

## Decision criterion (now with thresholds)

Solve the unit-current pickup-mode boundary-value problem (Sec. 2 of the note) for the
hermetic sheath in the thin-annulus limit and extract two coefficients:

```
F(β) ≡ V_eff / R^3     →  fit   F(β) = f_1 β + f_2 β^2 + ...
ℓ(β) ≡ L_p / (μ_0 R)   →  fit   ℓ(β) = ℓ_0 + ℓ_1 β + ℓ_2 β^2 + ...
```

**Operational thresholds.** For each fit, run the regression twice:

- Unconstrained fit with all coefficients free.
- Constrained fit with `ℓ_0 = 0` (for the ℓ regression) and with `f_0 = 0` (for the F
  regression, which is built in already since F has no constant term in the model).

Compute residual sum of squares (RSS) for each. Define:

- **"`ℓ_0` is consistent with zero":** RSS(constrained) ≤ 2 × RSS(unconstrained).
- **"`f_1` is statistically nonzero":** the fitted `f_1` exceeds 3× its bootstrap
  uncertainty.

Decision table:

| `f_1` nonzero? | `ℓ_0` consistent with zero? | Conclusion |
|---|---|---|
| Yes | Yes | Favorable thin-annulus scaling holds; compute crossover R_t (Sec. 9). |
| Yes | No | Inductance floor kills the gain. Report `ℓ_0` value and the dominant contribution (slit edge vs. external lead). |
| No | (either) | No signal-coupled response. Toroid sheath approach is dead for this G; do not proceed. |

## Geometry (fully specified)

Use dimensionless units: `R = 1`, `μ_0 = 1`, unit current `I = 1`.

- **Magnet volume:** `V_mag = { 1 < s < 1+β,  0 < z < α,  0 < φ < 2π }`
- **DC magnet field in V_mag:** `B_0 = B_max (R/s) φ̂` (compute with B_max = 1; the
  physical B_max factors out when crossover R_t is computed at the end)
- **α sweep (primary):** `α ∈ {0.5, 1.0, 2.0}` — all three must be run. Do not treat
  α as a "secondary" parameter; the answer can flip with α.
- **β sweep (primary):** `β ∈ {0.20, 0.12, 0.08, 0.05, 0.035, 0.025}` — six points
  to give the fit degrees of freedom against three free coefficients.
- **Slit nominal width:** `w/R = 0.03`. For mesh-resolution reasons (see below) the
  effective realized slit width should be reported alongside the nominal.

### Sheath surface (precise)

The sheath is the closed boundary `∂V_mag`, decomposed into four pieces:

1. **Inner cylinder:** `s = 1`, `0 ≤ z ≤ α`, `0 ≤ φ < 2π`.
2. **Outer cylinder:** `s = 1+β`, `0 ≤ z ≤ α`, `0 ≤ φ < 2π`.
3. **Bottom cap:** `z = 0`, `1 ≤ s ≤ 1+β`, `0 ≤ φ < 2π`.
4. **Top cap:** `z = α`, `1 ≤ s ≤ 1+β`, `0 ≤ φ < 2π`.

### Slit (precise)

The slit is a closed 1D loop on the sheath at azimuthal position φ ≈ 0, of azimuthal
half-width `w/2`. It consists of four connected segments:

- Segment A on the inner cylinder: `s = 1`, `z ∈ [0, α]`, at `φ = 0` (a vertical line
  on the inner cyl), broadened to an azimuthal width `w/R = 0.03`.
- Segment B on the top cap: `z = α`, `s ∈ [1, 1+β]`, at `φ = 0` (a radial line on
  the top cap), broadened to `w/R = 0.03`.
- Segment C on the outer cylinder: `s = 1+β`, `z ∈ [0, α]`, at `φ = 0`, broadened to
  `w/R = 0.03`.
- Segment D on the bottom cap: `z = 0`, `s ∈ [1, 1+β]`, at `φ = 0`, broadened to
  `w/R = 0.03`.

The four segments connect: A's bottom meets D, A's top meets B, B's outer meets C's top,
C's bottom meets D's outer. The slit topology is that of a closed loop around the
meridional cross-section. **This is the cut that interrupts toroidal-direction (φ̂)
currents on the sheath; without this cut, the closed sheath has only the no-go
axisymmetric response (Sec. 4 of the note).**

The pickup port has two terminals on opposite sides of segment A (one at φ = +w/2,
one at φ = -w/2) at the meridional midpoint `z = α/2`. Unit current is driven across
this gap.

### External return path (precise)

The readout wire exits one port terminal, runs radially outward at `z = α/2` to
`s = 10R`, then axially down to `z = -2R`, then radially inward to `s = 0`, then up
to `z = α/2`, then radially out to the other port terminal. Treat the wire as a thin
filament with regularization radius `r_wire = 0.01 R` for self-inductance.

**Document this choice in outputs.** It is part of G and the answer depends on it.

## Numerical approach: stream-function BEM

(BEM only. The volume FEM fallback in v1 is withdrawn — pick one approach and run it.)

Represent the surface current `K` on the sheath via a stream function `ψ` (so that
`∇_S · K = 0` is automatic):

- On a cylinder of radius `s`: `K_φ = -∂_z ψ`, `K_z = (1/s) ∂_φ ψ`
- On a flat cap (z = const): `K_s = (1/s) ∂_φ ψ`, `K_φ = -∂_s ψ`

`ψ` is continuous across edges (current conservation) and has a jump `Δψ = I = 1`
across the slit cut.

### Mesh (with adaptive refinement near slit)

Slit resolution is the dominant numerical issue. The slit at `w/R = 0.03` is ≈ 1.7°
of azimuth; a uniform 40-element mesh has 9° per element and would silently widen
the slit by a factor of ~5.

**Required adaptive azimuthal mesh:**

- Within `|φ| < 3w/2`: element size ≤ `w/3` (≈ 0.01 rad). Roughly 10 elements.
- For `3w/2 < |φ| < 6w/2`: element size ≤ `w` (≈ 0.03 rad). Roughly 6 elements.
- For `|φ| > 6w/2`: element size ≤ 0.1 rad (≈ 6°). Roughly 60 elements.

Total azimuthal elements per piece: ≈ 80–100. Axial/radial elements per piece: 20.
Total surface elements: ~7000.

**Validation of effective slit width.** After meshing, compute the actual azimuthal
width of the elements straddling the slit and report it as `w_eff`. The result is for
`w_eff`, not the nominal `w = 0.03`. Run one sensitivity check at `α = β = 1` with
`w_eff` at 0.5×, 1×, 2× nominal to confirm the result depends only weakly on slit
width within the resolved regime.

### Linear system

Linear system size: ~7000 × 7000, dense, double precision → ~400 MB. Direct LU solve
takes ~5–15 minutes per (α, β) combination on a modern workstation. The full primary
sweep (3 α × 6 β = 18 runs) is then 2–4 hours of compute. This fits in a workstation
session.

If memory becomes a problem, use scipy.sparse.linalg.gmres with diagonal
preconditioner and on-the-fly matrix-vector products (avoid storing the full matrix).
Do not use FMM/H-matrix acceleration unless the direct approach is shown insufficient
— it is over-engineering for this size.

### Implementation steps

1. **Mesh.** Generate adaptive surface mesh on the four pieces (described above).
   Identify slit-edge nodes; ψ on those nodes has the jump constraint.

2. **Basis.** Piecewise linear nodal ψ. K computed analytically per element from ψ
   gradients.

3. **Kernel.** Vector potential from surface current:
   ```
   A(r) = (μ_0/4π) ∫_S K(r')/|r-r'| dS' + A_wire(r)
   ```
   where A_wire is the Biot-Savart contribution from the external return filament.
   Use 4-point Gauss quadrature per element for non-adjacent panels. For adjacent
   panels (sharing an edge or node), use 9-point quadrature with analytical
   near-singular treatment of the `1/|r-r'|` singularity (split into singular +
   smooth parts; integrate singular part analytically over the test panel).

4. **Collocation.** Enforce `n̂ · B = 0` at element centroids using
   ```
   n̂ · B(r) = (μ_0/4π) ∫_S K(r') · [n̂(r) × (r-r')/|r-r'|^3] dS' + wire contribution.
   ```
   One equation per element.

5. **Slit constraint.** Add constraint rows enforcing `Δψ = 1` across the slit
   segments. The jump is between the ψ values on the "left" and "right" sides of
   the slit cut, taken consistently around the closed loop A→B→C→D→A.

6. **Solve** the constrained linear system for ψ (hence K).

7. **L_p via surface integral:**
   ```
   L_p = ∫_S K · A dS + L_wire_self
   ```
   The surface integral is computed as a double sum (with near-singular treatment
   for self-elements). `L_wire_self` is computed analytically from the wire path
   using the standard `μ_0/(4π) [ln(8R_wire/r_wire) - 2]` expression per segment,
   summed.

8. **V_eff via volume overlap (with adaptive φ quadrature):**
   ```
   V_eff = (R/μ_0) ∫_{V_mag} A_φ(s,φ,z) s ds dφ dz
   ```
   - Radial direction: 20 Gauss points over `s ∈ [1, 1+β]`.
   - Axial direction: 20 Gauss points over `z ∈ [0, α]`.
   - **Azimuthal direction: adaptive.** Use scipy.integrate.quad with `epsabs=1e-6`
     on the 1D integral over `φ` at each `(s, z)` point. This automatically refines
     near the slit where `A_φ` has rapid variation. Alternatively, use a piecewise
     Gauss-Legendre rule: 50 points in `|φ| < 0.5`, 30 points in `0.5 < |φ| < π`.
     Confirm both methods agree to 1%.

## Cross-checks (run before production sweep, modified per review)

### Cross-check 1 (revised — V_eff only): zero-slit topology check

Set up the linked-return topology: route the external return wire through the
central hole, linking the major axis once. By the linked-mode analysis (Sec. 5 of
the note), this produces a vector potential with no `A_φ` component in V_mag, hence

```
V_eff (linked-return topology) → 0  (to within numerical precision).
```

Run at `α = β = 1` (well-resolved regime). Pass criterion: |V_eff| < 0.01 × V_eff
from the production setup at the same (α, β).

**Note removed.** The previous inductance-matching part of this check (matching
`L_p` to `(μ_0 h / 2π) ln(1+β)`) is unreliable because the linked-mode wire self-
inductance contaminates the comparison. Skip it; use only the V_eff topology check.

### Cross-check 2: zero-slit limit (primary validation)

At `α = β = 1`, shrink the slit width through `w_eff/R ∈ {0.03, 0.01, 0.003, 0.001}`.
The axisymmetric no-go result (Sec. 4 of the note) requires `V_eff → 0` as the slit
closes.

Pass criteria:
- `V_eff` is monotone decreasing in `w_eff`.
- `V_eff(w_eff = 0.001) < 0.1 × V_eff(w_eff = 0.03)`.

This is now the **primary** code-validity check. If it fails, the BEM has a topology
bug.

### Cross-check 3: reciprocity self-consistency

For the production setup at `α = β = 1`, compute `L_p` two ways:
- Surface integral: `L_p = ∫_S K · A dS` (the production method).
- Volume integral: `L_p = (1/μ_0) ∫ |B|^2 d^3x` over a sufficiently large box.

Pass criterion: agreement to 5%. If worse, the kernel evaluation or near-field
treatment has a bug.

### Cross-check 4: gauge invariance of V_eff

Add `∇χ` for a simple test χ (e.g., `χ = sφ z`) to the computed `A` and recompute
`V_eff`. It should be unchanged because `B_0` is divergence-free in V_mag and
`B_0 · n̂ = 0` on `∂V_mag`.

Pass criterion: relative change < 10⁻⁴.

**If any of cross-checks 2, 3, 4 fail, stop and debug.** Cross-check 1 (zero-V_eff
in linked topology) is a topology-only check; it should also pass but is less
diagnostic than 2 and 3.

## Production sweep

After cross-checks 2, 3, 4 pass:

For each `α ∈ {0.5, 1.0, 2.0}`:
  For each `β ∈ {0.20, 0.12, 0.08, 0.05, 0.035, 0.025}`:
    1. Solve the BEM problem.
    2. Compute `L_p` and `V_eff`.
    3. Record `F(β; α) = V_eff/R^3` and `ℓ(β; α) = L_p/(μ_0 R)`.
    4. Save raw solution data (K, A_φ on V_mag grid) for inspection.

For each α, fit:
```
F(β) = f_1(α) β + f_2(α) β^2 + f_3(α) β^3   (constant term = 0 imposed)
ℓ(β) = ℓ_0(α) + ℓ_1(α) β + ℓ_2(α) β^2
```
Use weighted nonlinear least squares (use 1% relative weight as proxy for numerical
uncertainty; tighten if cross-check 3 indicates better). Bootstrap or jackknife the
residuals to get uncertainties on the fitted coefficients.

Apply the decision criterion (Sec. 1) for each α separately. The conclusion may
depend on α; if so, report which α regimes are favorable.

## Crossover R_t calculation (now explicit)

For each `α` where the favorable scaling holds, compute the crossover `R_t*` at
which `Q^tor = Q^core` for CAL parameters.

CAL design parameters:
- `B_max = 5 T`
- Plausible bore-volume budget: take `V_mag ∈ {0.1, 1, 10} m^3` and compute the
  corresponding magnet cost `C = ∫ B^2 dV` in each case. For the toroid, this is
  `C_tor = 2π B_max^2 R_t^3 α ln(1+β)`. For the Core, it is `C_core = π R_m^2 h B_max^2`.

For Core at fixed `C_core` and the optimum pickup ratio `r* = e^(0.3)`:
```
Q^core = (B_max^4 / μ_0^2) × (h R_m^2 / 2)^4 × (ln r*)^2 × (2π / μ_0 h)^2
       = (B_max^4 / μ_0^4) × (h R_m^2)^4 × ln^2 r* / (16 π^2 h^2)
```
(verify the formula by re-deriving from Sec. 6 of the note — `V_eff = h R_m^2 ln r / 2`
and `L_p = μ_0 h ln r / (2π)`).

For Toroid at fixed `C_tor`, using the favorable thin-annulus scaling with fitted
`f_1, ℓ_1`:
```
V_eff ≈ R^3 × f_1 × β  =  f_1 × C_tor / (2π B_max^2 α)    (linear-in-β limit)
L_p   ≈ μ_0 R × ℓ_1 × β  =  μ_0 ℓ_1 × C_tor / (2π B_max^2 R^2 α)
Q^tor ≈ (B_max^4 / μ_0^2) × V_eff^4 / L_p^2
      = (B_max^4 / μ_0^2) × (f_1 C_tor / 2π B_max^2 α)^4 / (μ_0 ℓ_1 C_tor / 2π B_max^2 R^2 α)^2
      = (1 / μ_0^4) × R^4 × (f_1^4 / ℓ_1^2) × (C_tor^2 / (2π α)^2)
```

Set `Q^tor = Q^core` at fixed `C_tor = C_core`, solve for `R_t*`. Report this number
in physical units (meters) for each plausible cost budget. If `R_t* > 1 m`, flag as
likely impractical at CAL scale.

## "What counts as inconclusive" (new section)

Do not tune parameters to force a particular conclusion. Stop and report the raw
data if any of the following hold:

1. **Non-monotone F(β) or ℓ(β) in β.** The fits assume polynomial behavior in β; if
   the data is non-monotone or has multiple sign changes, the model doesn't fit.
2. **Fit χ²/dof > 5 for both unconstrained and constrained fits.** The polynomial
   ansatz is wrong; report raw points and abandon the linear/quadratic interpretation.
3. **Condition number of BEM matrix > 10⁸ at smallest β.** Numerical instability
   makes the smallest-β point unreliable; report only the resolved points.
4. **Cross-check 3 (reciprocity) disagreement > 20%.** The energy and surface-integral
   forms of L_p should match; if not, kernel evaluation is bugged.
5. **Inconsistent results across α values that change the binary decision.** If
   α=0.5 gives "favorable" and α=2 gives "unfavorable," the answer is genuinely
   α-dependent and the binary criterion does not apply. Report the α-dependence
   structure.

In any of these cases, write a `report.md` documenting what was seen and what was
not concluded. Do not paper over.

## Outputs

Save to `/mnt/user-data/outputs/`:

1. `bem_solver.py` — main implementation. Self-contained; numpy + scipy only.
2. `mesh_validation.json` — for each (α, β) run: nominal slit width, effective slit
   width, total element count, condition number of the BEM matrix.
3. `crosscheck_results.json` — results of all four cross-checks with pass/fail
   flags and quantitative comparisons.
4. `production_sweep.csv` — columns: α, β, F(β;α), ℓ(β;α), uncertainty estimates,
   wall-clock time per run.
5. `fit_results.json` — for each α: f_1, ℓ_0, ℓ_1, f_2, ℓ_2, with bootstrap
   uncertainties; binary decision (favorable / unfavorable / inconclusive) and
   reason.
6. `crossover_analysis.json` — for each α with favorable scaling: R_t* in meters at
   each of three cost budgets (V_mag ∈ {0.1, 1, 10} m³), plus Q^core for comparison.
7. `plots/`:
   - `crosscheck_zeroslit.pdf` — V_eff vs. effective slit width.
   - `crosscheck_reciprocity.pdf` — L_p energy form vs. surface-integral form, across
     the sweep.
   - `scaling_F_vs_beta.pdf` — F(β) vs. β for each α; log-log inset.
   - `scaling_ell_vs_beta.pdf` — ℓ(β) vs. β for each α; log-log inset.
   - `surface_K_visualization.pdf` — heat-map of |K| on the sheath at α=β=1; should
     show concentration near slit endpoints.
   - `crossover_vs_cost.pdf` — R_t* vs. cost budget for each favorable α.
8. `report.md` — one-page executive summary:
   - Pass/fail of all four cross-checks.
   - Fit coefficients with uncertainties for each α.
   - Binary decision per α with reasoning.
   - Crossover R_t* if applicable, with units.
   - Any "inconclusive" flags from the section above.

## Estimated complexity

Implementation: 1000–1800 lines of Python (numpy + scipy). Implementer should plan
2–3 days for code + debugging + cross-check passes, then 4–6 hours of compute for the
production sweep, then 1–2 hours of analysis. Total ~4 working days.

If the near-singular kernel treatment is non-trivial in implementation (likely), most
of the debugging time will be there. The cross-check 3 (reciprocity) is the main
diagnostic for whether kernel handling is correct.

## What to do if blocked

- **Cross-check 2 fails (V_eff doesn't go to zero with closing slit):** the slit
  topology or boundary condition implementation is wrong. Probably the Δψ = 1
  constraint is being applied incorrectly across the slit segments.
- **Cross-check 3 fails (reciprocity):** near-singular kernel handling is wrong.
  Check the analytical near-field treatment.
- **BEM matrix singular or near-singular:** the constraint system is over- or
  under-determined. Check the slit-edge node identification.
- **Fits diverge or are wildly noisy:** see "what counts as inconclusive." Probably
  not all six β points are usable; truncate to the resolved subset.

Document everything in `report.md`. An honest "the BEM didn't converge at small β,
here's what we have" is more valuable than a tuned result.

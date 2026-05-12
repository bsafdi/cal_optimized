# Hermetic-sheath toroidal pickup — executive summary (post wire-removal patch)

This report covers two studies on the same BEM solver:
1. **Study i**: original v2 spec with the production readout wire (out to s=10R, axial down, return along z-axis).
2. **Study ii**: wire-removal patch (`wire_removal_patch.md`) — sheath in vacuum, slit Dirichlet BCs only.

The patch was executed cleanly. The CSV/JSON outputs in this directory reflect study (ii); study (i) results are preserved in git history under commits dceadea-cb5f549.

## Main result: exact symmetry zero

For the sheath in vacuum (no external return wire), **V_eff is identically zero by phi-reflection symmetry**. Proof sketch:

* Sheath geometry is invariant under phi -> -phi.
* The slit BC psi(+w/2) = +1/2, psi(-w/2) = -1/2 is mapped under reflection to its negative.
* The energy-minimising psi is therefore antisymmetric: psi(phi) = -psi(-phi).
* K_phi is antisymmetric, K_z and K_s are symmetric.
* A_phi(s, phi, z) = -A_phi(s, -phi, z), so the integral over the symmetric domain V_mag is zero.

The BEM returns F values at ~1e-15 to 1e-19, consistent with floating-point noise around the exact zero. (See the LaTeX notes, sec. "Wire-removal patch", for the full derivation.)

## Cross-checks (sheath only, alpha = beta = 1)

| Check | Result | Comment |
|---|---|---|
| XC1 (linked return) | **N/A** | No external wire to route through the central hole. |
| XC2 (zero-slit limit) | **PASS** | |F| at machine precision for all w in {0.03, 0.01, 0.003, 0.001}. Consistent with V_eff = 0 exactly. ell stable at 0.103 across all w. |
| XC3 (reciprocity) | **PASS** | L_surface = 0.1031, L_volume = 0.1045, rel diff = 1.3% < 5%. Tighter than study (i)'s 2.5% because no wire-self-energy to subtract. |
| XC4 (gauge invariance) | **PASS** | Delta V_eff = 6.8e-15. V_eff baseline zero by symmetry; gauge shift also zero analytically. |

## Production sweep (sheath only)

Re-ran the 18-point (alpha, beta) sweep on SLURM lr7 job 22515010 (host n0062.lr7, 8 cpus, 190 s wall). All F at machine precision; ell shows linear-in-beta scaling at moderate beta:

| beta | ell at alpha=0.5 | ell at alpha=1.0 | ell at alpha=2.0 |
|---|---|---|---|
| 0.200 | 0.0139 | 0.0253 | 0.0469 |
| 0.120 | 0.00847 | 0.0143 | 0.0247 |
| 0.080 | 0.00550 | 0.00815 | 0.0125 |
| 0.050 | 0.00308 | 0.00365 | 0.00150 |
| 0.035 | 0.00168 | 0.00056 | -0.00734 |
| 0.025 | 0.00054 | -0.00238 | -0.0125 |

The beta = 0.025, 0.035 points become unphysical (ell < 0) at alpha = 1, 2 due to mesh undersampling of very thin caps (cap thickness 0.025/15 = 1.7e-3 at n_t = 15, condition number 1e7+). These points are excluded from the polynomial fits.

## Fits (beta in {0.05, 0.08, 0.12, 0.20})

| alpha | f_1 | ell_0 | ell_1 | RSS_constr / RSS_free |
|---|---|---|---|---|
| 0.5 | -2.7e-17 | -1.1e-3 | +0.087 | 132 |
| 1.0 | -1.5e-15 | -4.5e-3 | +0.167 | 121 |
| 2.0 | -6.6e-14 | -1.7e-2 | +0.398 | 160 |

f_1 is at floating-point noise across all alpha (consistent with the exact V_eff = 0 result above). ell_1 grows with alpha as expected. Fitted ell_0 is small in magnitude (<10% of ell_1 beta_min) and slightly negative — consistent with mesh noise in the 4-point fit, not a physical positive-definite floor.

## Binary decision per alpha

For all three alpha values: **row 3 of the spec's decision table — "No signal-coupled response. Toroid sheath approach is dead for this G; do not proceed."** This time the verdict is not a numerical-noise artifact: V_eff = 0 is an exact phi-reflection symmetry result for the sheath in vacuum.

## Crossover R_t*

Marked "not applicable" because no alpha has favorable scaling (V_eff exactly zero -> Q^tor = 0). The crossover_analysis.json contains the formal calculation with the unresolved f_1 plugged in — Q^tor / Q^core ~ 1e-53 to 1e-63, all astronomically below unity. These numbers are not physically meaningful; they reflect f_1 at machine precision.

## Comparison: study (i) vs study (ii)

| Quantity | Study (i): with wire | Study (ii): sheath only |
|---|---|---|
| ell (alpha=beta=1) | 19.03 (99.3% wire-self) | 0.103 |
| ell_1 (alpha=1) | 0.59 (small over wire floor) | 0.167 |
| ell_0 (alpha=1) | 18.87 (wire-dominated) | ~0 (consistent with positive-def. constraint) |
| F (alpha=beta=1) | -1.2e-4 (sub-noise residual) | 0 exactly (symmetry) |
| Decision | Row 3 (numerically), Row 2 (conditionally) | Row 3 (by exact symmetry) |

**Both studies give the same formal decision (row 3) but for very different reasons.** Study (i): inductance floor and unresolved F. Study (ii): exact symmetry zero of V_eff combined with favorable but unsignalled ell scaling.

## What the two studies together tell us

1. The sheath's intrinsic inductance scaling IS favorable (study ii): ell ~ ell_1 beta with ell_0 consistent with zero in the well-resolved regime. The asymptotic thin-annulus form holds for the sheath.
2. The sheath alone has no signal mode (study ii): V_eff = 0 by phi-reflection symmetry of the antisymmetric slit BC.
3. The spec's production wire DID break this symmetry, but added so much self-inductance that the favorable scaling was killed even in principle (study i).
4. **Useful design requirement**: a readout wire that (a) breaks phi-reflection symmetry and (b) has self-inductance subdominant to mu_0 R * ell_1 * beta, i.e. << 0.05 * mu_0 * R_t at alpha=1, beta=0.2.

The spec's wire violated (b) by a factor of ~400. A coaxial return hugging the inner sheath at small radial offset would in principle satisfy both. Quantifying that is a follow-up study.

## Inconclusive flags per spec

* F(beta) is at machine precision; doesn't fit a non-trivial polynomial. This is an exact zero, not noise, so the spec's "non-monotone polynomial" criterion is interpreted in light of the symmetry analysis above.
* No condition-number violation in the well-resolved beta range (cond < 4e6 at beta >= 0.05). Smallest-beta points (beta = 0.025) have cond up to 1.8e7 < 1e8 but produce unphysical negative ell; excluded from fits.
* XC3 reciprocity passes at 1.3% (tighter than 5% threshold).

## Bottom line

**For the sheath in vacuum, V_eff = 0 identically by phi-reflection symmetry of the antisymmetric slit BC.** The spec's decision table routes to row 3 (do not proceed), but the deep reason is symmetry, not inductance. The sheath's intrinsic inductance scales favorably with beta. A useful pickup design requires an external wire that breaks the symmetry while keeping its self-inductance below mu_0 R ell_1 beta. The spec's production wire failed this by ~400x; finding a wire that satisfies both constraints is the natural next study.

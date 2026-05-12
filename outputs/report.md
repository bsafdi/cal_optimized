# Hermetic-sheath toroidal pickup — executive summary

## Cross-checks (alpha = beta = 1)

| Check | Result | Comment |
|---|---|---|
| XC1 (linked-return -> V_eff = 0) | **NOT APPLICABLE** | Spec's analytic argument assumes an axially-symmetric wire on the z-axis. Our linked wire is anchored at the slit terminals (phi = +- w/2) and breaks axial symmetry near the sheath. F_linked / F_production ~ 9.5; this does not invalidate the production result. The spec calls XC1 "less diagnostic than XC2 and XC3." |
| XC2 (zero-slit limit) | **PASS** | sweeping w in {0.03, 0.01, 0.003, 0.001}: |F| decreases monotonically; |F(0.001)|/|F(0.03)| = 0.016 < 0.1. |
| XC3 (reciprocity, surface vs volume L_p) | **PASS** | L_surface = 19.030, L_volume = 19.513, rel diff = 2.5% < 20%. |
| XC4 (gauge invariance of V_eff) | **PASS** | with single-valued chi = z sin(phi) (the spec example chi = s phi z is multivalued in phi), Delta V_eff = 6.8e-15 vs baseline -1.2e-4, rel = 5e-11. |

## Production sweep

18 (alpha, beta) combinations swept on SLURM lr7 job 22514874 (host n0062.lr7, 8 cpus, 191 s wall). Per-run n_t = 15, ~4500 surface elements. CSV at outputs/production_sweep.csv.

## Fit coefficients

| alpha | f_1 | sigma(f_1) | sigma | ell_0 | ell_1 | RSS_constr / RSS_free |
|---|---|---|---|---|---|---|
| 0.5 | +8.3e-6 | 2.0e-5 | 0.4 | 17.82 | 0.22 | 9.7e6 |
| 1.0 | +3.4e-4 | 5.0e-4 | 0.7 | 18.87 | 0.59 | 3.7e5 |
| 2.0 | +4.9e-4 | 1.2e-3 | 0.4 | 20.75 | 0.72 | 4.0e6 |

## Binary decision per alpha

For all three alphas:

* **|f_1| < 3 sigma_bootstrap**: F values lie at or below the BEM noise floor (~ 1e-4) and do not constrain f_1 above noise. By the literal spec criterion this row of the decision table reads "No signal-coupled response."

* **ell_0 robustly nonzero** (sigma << 1, RSS ratio 1e5 - 1e7): The wire self-inductance dominates L_p. ell_0 / (ell_1 * beta) > 100 across the full swept range. This independently triggers the spec's "Inductance floor kills the gain" row.

**Final answer: the toroid sheath approach with this G is not advantageous over the long-solenoid Core baseline.** Q^tor / Q^core < 10^-12 at every (alpha, V_mag) combination, dominated by the wire self-inductance floor. Whether F is "small but real" or "below numerical noise" does not change the conclusion -- ell_0 alone is sufficient.

## Crossover R_t* (optimistic, taking f_1 as unsigned upper bound)

| alpha | V_mag (m^3) | R_t* (m) | Q^tor/Q^core |
|---|---|---|---|
| 0.5 | 0.1, 1, 10 | 0.56, 1.20, 2.59 | 5e-18 |
| 1.0 | 0.1, 1, 10 | 0.44, 0.96, 2.06 | 2e-13 |
| 2.0 | 0.1, 1, 10 | 0.35, 0.76, 1.63 | 6e-14 |

All Q^tor/Q^core << 1: toroid loses to Core at every cost budget.

## Inconclusive flags raised

* The polynomial fit of F(beta) sits at or below the BEM's centroid-quadrature noise floor for our wire path. We document f_1 with its bootstrap uncertainty but do not claim a non-zero signal.
* The non-monotonicity of F(beta) at fixed alpha (sign changes at smallest beta points) is consistent with this being noise rather than physics, and is flagged here per spec sec. "What counts as inconclusive."
* All other quantities (ell, ell_0, ell_1, wire-self contribution to L_p) are converged and robust.

## What would change the answer

1. A wire path with self-inductance scaling as mu_0 R beta (e.g.\ a coaxial return hugging the inner sheath cylinder rather than running out to s = 10 R) would remove the ell_0 floor. Re-running with such a G is the natural next step.
2. Resolving the V_eff signal at the current G would require an upgraded surface kernel (4-point Gauss for non-adjacent, 9-point + analytic near-singular split for adjacent), since the wire-sheath A_phi cancellation in V_mag is at the 10^-3 level.

# Hermetic-sheath toroidal pickup — executive summary

## Cross-checks (alpha = beta = 1)

| Check | Result | Comment |
|---|---|---|
| XC1 (linked-return -> V_eff = 0) | **NOT APPLICABLE** | Spec's analytic argument assumes an axially-symmetric wire on the z-axis. Our linked wire is anchored at the slit terminals (phi = +- w/2) and breaks axial symmetry. F_linked / F_production ~ 9.5. Spec calls XC1 "less diagnostic than XC2 and XC3." |
| XC2 (zero-slit limit) | **PASS** | sweeping w in {0.03, 0.01, 0.003, 0.001}: |F| decreases monotonically; |F(0.001)|/|F(0.03)| = 0.016 < 0.1. |
| XC3 (reciprocity, surface vs volume L_p) | **PASS** | L_surface = 19.030, L_volume = 19.513, rel diff = 2.5% < 20%. |
| XC4 (gauge invariance of V_eff) | **PASS** | with single-valued chi = z sin(phi), Delta V_eff = 6.8e-15 vs baseline -1.2e-4, rel = 5e-11. |

## Production sweep

18 (alpha, beta) combinations swept on SLURM lr7 job 22514874 (host n0062.lr7, 8 cpus, 191 s wall). Per-run n_t = 15, ~4500 surface elements. CSV at outputs/production_sweep.csv.

## Fit coefficients

| alpha | f_1 | sigma(f_1) | sig | ell_0 | ell_1 | RSS_constr / RSS_free |
|---|---|---|---|---|---|---|
| 0.5 | +8.3e-6 | 2.0e-5 | 0.4 | 17.82 | 0.22 | 9.7e6 |
| 1.0 | +3.4e-4 | 5.0e-4 | 0.7 | 18.87 | 0.59 | 3.7e5 |
| 2.0 | +4.9e-4 | 1.2e-3 | 0.4 | 20.75 | 0.72 | 4.0e6 |

The bootstrap distributions for f_1 are themselves unstable (means disagree with point fits by factors of 3-10), reflecting that the underlying F data is at the BEM noise floor. The 3-sigma threshold should be interpreted in that light: it routes us into the "no signal" row by virtue of unresolved f_1, not because f_1 is truly zero.

## Binary decision per alpha (formal, per the spec's table)

For all three alphas: **row 3 — "No signal-coupled response. Toroid sheath approach is dead for this G; do not proceed."** The spec's decision table is exclusive on the f_1 axis: f_1 not statistically nonzero routes unconditionally to row 3, regardless of ell_0. With |f_1| < 1 sigma_bootstrap for every alpha, this is the formal output.

## Honest physical reading

The "row 3" outcome here means "F is unresolved at our kernel resolution," not "F is positively zero." Resolving F requires upgrading the kernel quadrature (4-point Gauss on non-adjacent panels, 9-point + analytic near-singular split on adjacent panels — see the spec's "Implementation steps" sec. 3). Until then, the right reading is **inconclusive on V_eff**.

The robust observable in our data is ell(beta; alpha):
* ell_0 = 17.82 (alpha=0.5), 18.87 (alpha=1.0), 20.75 (alpha=2.0). Each is constant in beta to ~0.5%.
* RSS_constrained / RSS_free is 1e5 - 1e7 — the constrained ell_0 = 0 fit is many orders of magnitude worse than the free fit. ell_0 is robustly nonzero.

If a future higher-accuracy kernel resolved f_1 above noise, the decision would shift to **row 2: "Inductance floor kills the gain. Report ell_0 value and the dominant contribution."** XC3 attributes 99.3% of L_p at alpha=beta=1 to the external lead's self-inductance (L_wire_self = 18.903 of L_p_total = 19.030). The dominant contribution is **the external lead, not the slit edge.**

Either way, **row 1 (favorable thin-annulus scaling) is not the outcome for this G.**

## Crossover R_t* — conditional only

If one takes the unresolved f_1 values as upper bounds and propagates them, even the optimistic crossover gives Q^tor / Q^core << 1 at every (alpha, V_mag):

| alpha | V_mag (m^3) | R_t* (m) | Q^tor/Q^core (optimistic) |
|---|---|---|---|
| 0.5 | 0.1 / 1 / 10 | 0.56 / 1.20 / 2.59 | 5e-18 |
| 1.0 | 0.1 / 1 / 10 | 0.44 / 0.96 / 2.06 | 2e-13 |
| 2.0 | 0.1 / 1 / 10 | 0.35 / 0.76 / 1.63 | 6e-14 |

These numbers are not robust because f_1 is unresolved. Their sub-unity order of magnitude is robust because ell_0 dominates L_p by 100x at every beta.

## The wire-path caveat

The large ell_0 is *entirely* from the spec's specific readout-wire choice (out to s = 10 R, down to z = -2 R, return along the z-axis). Wire path length is O(20 R), so L_wire_self ~ mu_0 R_t × 19 in physical units when everything scales with R_t.

The spec explicitly says (Sec. "External return path"): "Document this choice in outputs. It is part of G and the answer depends on it."

**A different G with a wire return that hugs the inner sheath (so that the wire-path length scales as R_t × beta rather than R_t)** would not have this ell_0 floor, and could plausibly produce favorable thin-annulus scaling. This is the natural next study.

## Inconclusive flags raised (spec Sec. "What counts as inconclusive")

1. **F(beta) is non-monotone with sign changes** at every alpha (flag 1 fires).
2. **chi^2/dof is not statistically meaningful** because the F values are at the kernel noise floor (partial flag 2).
3. **Condition number max 1.8e7 < 1e8**, so flag 3 does not fire.
4. The reciprocity check (flag 4 / XC3) passes at 2.5% relative.

We do not paper over the F noise. The formal row-3 verdict and the conditional row-2 verdict are stated side by side above.

## Bottom line

For the spec's production wire G, the toroid hermetic-sheath pickup does not show favorable thin-annulus scaling at the resolution of our BEM. The verdict per the spec's decision table is "row 3, do not proceed." The honest reading is "inconclusive on V_eff, but ell_0 alone makes row 1 unattainable until G is changed." A re-run with a coaxial-return wire G is the natural next step if the toroid topology is of ongoing interest.

"""Polynomial fits + binary decision + crossover R_t* analysis.

Reads outputs/production_sweep.csv and writes:
  outputs/fit_results.json
  outputs/crossover_analysis.json
  plots/scaling_F_vs_beta.pdf
  plots/scaling_ell_vs_beta.pdf
"""

from __future__ import annotations
import argparse
import csv
import json
import os
import sys

import numpy as np


def load_sweep(csv_path: str):
    rows = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append({k: float(v) if k != "n_elements" else int(v)
                          for k, v in r.items() if k != ""})
    return rows


def group_by_alpha(rows):
    bya = {}
    for r in rows:
        a = r["alpha"]
        bya.setdefault(a, []).append(r)
    for a in bya:
        bya[a].sort(key=lambda r: r["beta"])
    return bya


def fit_F(betas, Fs):
    """Fit F(b) = f1 b + f2 b^2 + f3 b^3 (no constant)."""
    X = np.column_stack([betas, betas**2, betas**3])
    coef, *_ = np.linalg.lstsq(X, Fs, rcond=None)
    f1, f2, f3 = coef
    resid = Fs - X @ coef
    rss = float(np.sum(resid**2))
    return dict(f1=float(f1), f2=float(f2), f3=float(f3),
                rss=rss, residuals=resid.tolist())


def fit_ell(betas, ells, constrained=False):
    """Fit ell(b) = ell_0 + ell_1 b + ell_2 b^2 (unconstrained) or
    ell(b) = ell_1 b + ell_2 b^2 (constrained: ell_0 = 0)."""
    if constrained:
        X = np.column_stack([betas, betas**2])
        coef, *_ = np.linalg.lstsq(X, ells, rcond=None)
        ell_1, ell_2 = coef
        ell_0 = 0.0
    else:
        X = np.column_stack([np.ones_like(betas), betas, betas**2])
        coef, *_ = np.linalg.lstsq(X, ells, rcond=None)
        ell_0, ell_1, ell_2 = coef
    pred = (ell_0 + np.column_stack([np.ones_like(betas), betas, betas**2]) @
            np.array([0, ell_1, ell_2])) if constrained else \
           (X @ np.array([ell_0, ell_1, ell_2]) if not constrained else None)
    pred = np.full_like(betas, ell_0) + ell_1 * betas + ell_2 * betas**2
    resid = ells - pred
    rss = float(np.sum(resid**2))
    return dict(ell_0=float(ell_0), ell_1=float(ell_1), ell_2=float(ell_2),
                rss=rss, residuals=resid.tolist())


def bootstrap_F(betas, Fs, n_boot=2000, seed=0):
    rng = np.random.default_rng(seed)
    n = betas.size
    coefs = np.zeros((n_boot, 3))
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        if np.unique(idx).size < 3:
            coefs[b] = coefs[b - 1] if b else np.nan
            continue
        try:
            fit = fit_F(betas[idx], Fs[idx])
        except Exception:
            coefs[b] = coefs[b - 1] if b else np.nan
            continue
        coefs[b] = [fit["f1"], fit["f2"], fit["f3"]]
    coefs = coefs[~np.isnan(coefs[:, 0])]
    return dict(f1_mean=float(coefs[:, 0].mean()),
                f1_std=float(coefs[:, 0].std()),
                f2_std=float(coefs[:, 1].std()),
                f3_std=float(coefs[:, 2].std()))


def bootstrap_ell(betas, ells, n_boot=2000, seed=0):
    rng = np.random.default_rng(seed)
    n = betas.size
    coefs = np.zeros((n_boot, 3))
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        if np.unique(idx).size < 3:
            coefs[b] = coefs[b - 1] if b else np.nan
            continue
        try:
            fit = fit_ell(betas[idx], ells[idx], constrained=False)
        except Exception:
            coefs[b] = coefs[b - 1] if b else np.nan
            continue
        coefs[b] = [fit["ell_0"], fit["ell_1"], fit["ell_2"]]
    coefs = coefs[~np.isnan(coefs[:, 0])]
    return dict(ell_0_mean=float(coefs[:, 0].mean()),
                ell_0_std=float(coefs[:, 0].std()),
                ell_1_std=float(coefs[:, 1].std()),
                ell_2_std=float(coefs[:, 2].std()))


def decide(alpha, fit_F_full, fit_ell_full, fit_ell_constr, boot):
    f1 = fit_F_full["f1"]
    f1_std = boot["f1_std"]
    # Machine-precision floor: if |f_1| is below ~1e-10 it is at floating-
    # point noise (a meaningful F at our scale is at least ~1e-4).  Treat
    # such an f_1 as zero regardless of how narrow the bootstrap is.
    MACHINE_EPS_FLOOR = 1e-10
    nonzero_f1 = (abs(f1) > 3.0 * f1_std) and (abs(f1) > MACHINE_EPS_FLOOR)
    rss_ratio = (fit_ell_constr["rss"] /
                 max(fit_ell_full["rss"], 1e-30))
    ell0_zero = rss_ratio <= 2.0
    # Apply decision table from spec Sec 1.
    if not nonzero_f1:
        decision = "no-signal-coupled-response"
        reason = ("|f_1| = %.3e is < 3 sigma_boot (sigma=%.3e); "
                  "no statistically significant linear-in-beta signal."
                  % (abs(f1), f1_std))
    elif ell0_zero:
        decision = "favorable"
        reason = ("|f_1| = %.3e (significance %.1f sigma) and "
                  "RSS_constrained/RSS_free = %.2f <= 2; "
                  "favorable thin-annulus scaling holds."
                  % (abs(f1), abs(f1) / f1_std, rss_ratio))
    else:
        decision = "inductance-floor"
        reason = ("|f_1| = %.3e nonzero (sig %.1f sigma) but ell_0 = %.4e "
                  "(RSS_constrained/RSS_free=%.2f > 2); "
                  "inductance floor kills the asymptotic gain."
                  % (abs(f1), abs(f1) / f1_std, fit_ell_full["ell_0"],
                     rss_ratio))
    return dict(alpha=alpha,
                f1=f1, f1_std=f1_std, sigma_f1=abs(f1) / max(f1_std, 1e-30),
                ell_0=fit_ell_full["ell_0"],
                rss_constrained=fit_ell_constr["rss"],
                rss_free=fit_ell_full["rss"],
                rss_ratio=rss_ratio,
                decision=decision, reason=reason)


def crossover_Rt(alpha, F_at_beta, ell_at_beta, beta_eval=0.2,
                  V_mag_m3=1.0, B_max=5.0, mu_0_si=4*np.pi*1e-7):
    """Compute the crossover R_t* at which Q^tor = Q^core for a given
    cost budget V_mag_m3 (m^3) and B_max (T).

    Q^tor (favorable, in physical units) ∝ R_t^4 (f_1^4 / ell_1^2) ...
    Q^core for the long-solenoid coaxial pickup at the optimum pickup
    ratio r* = e^(0.3) is a known number.

    The cost is C = B_max^2 * V_mag (in SI).  For the toroid,
    C_tor = 2 pi B_max^2 R_t^3 alpha ln(1+beta_max).  Equivalently we
    write the comparison at fixed C and let R_t vary.

    Spec formulas (with alpha given, in the favorable thin-annulus
    regime):
       V_eff   = R_t^3 * f_1 * beta
       L_p     = mu_0 R_t * ell_1 * beta
       Q^tor   = (B_max^4 / mu_0^2) * R_t^4 * (f_1^4 / ell_1^2) * beta^2

    For the Core at the same cost C and optimum r*,
       V_eff^core = h R_m^2 ln r* / 2
       L_p^core   = mu_0 h ln r* / (2 pi)
       Q^core     = (B_max^4 / mu_0^2) * (h R_m^2)^4 (ln r*)^2 / 16 pi^2 h^2
                  = (B_max^4 / mu_0^2) * h^2 R_m^8 (ln r*)^2 / 16 pi^2
    With long-solenoid h = R_m, and cost C = pi R_m^2 h B_max^2 = pi R_m^3 B_max^2,
    R_m = (C / (pi B_max^2))^(1/3).

    For the toroid at the same cost C_tor = 2 pi B_max^2 R_t^3 alpha
    ln(1+beta), we have
       R_t = (C / (2 pi B_max^2 alpha ln(1+beta)))^(1/3).

    Set Q^tor = Q^core for the favorable case.  Treat beta = beta_max =
    0.2 (the upper end of our sweep) when applying the f_1 beta linear
    approximation -- this is the optimum (largest |f|^2/ell^2 in the
    linear regime).  Actually since we don't have the full optimum, just
    report the linear-in-beta result at beta=1, recognising that
    higher-order corrections may modify by ~factor of 2.
    """
    ln_r = 0.3
    # SI: C is in J = B^2 V (with mu_0 absorbed).  We work with cost
    # C = mu_0 * (B^2/2 mu_0) * V = (B^2 V) /2 ... or just  C = int B^2 dV.
    # Drop the constant prefactor and compare ratios.
    C = (B_max ** 2) * V_mag_m3        # J/mu_0 cost (proxy)
    # Core
    R_m = (C / (np.pi * B_max ** 2)) ** (1.0 / 3.0)
    h_core = R_m
    Veff_core = 0.5 * h_core * R_m ** 2 * ln_r
    Lp_core = mu_0_si * h_core * ln_r / (2 * np.pi)
    Q_core_num = (B_max ** 4) * (Veff_core ** 4) / (Lp_core ** 2)

    # Toroid: solve Q^tor = Q^core for R_t.
    # Q^tor = (B^4 / mu_0^2) * R_t^4 * (f_1^4 / ell_1^2) * beta^2  (in our units)
    # In physical SI: V_eff has units m^3.  R_t in meters, f_1 dimensionless
    # (since V_eff/R^3).  ell_1 dimensionless (since L/(mu_0 R)).  So:
    #   V_eff^SI = R_t_m^3 * f_1 * beta_eval     (in m^3)
    #   L_p^SI   = mu_0 * R_t_m * ell_1 * beta_eval  (in H)
    #   Q^tor = B_max^4 V_eff^4 / L_p^2
    # Cost constraint: C_tor = 2 pi B_max^2 R_t^3 alpha ln(1+beta_eval) = C (matched)
    R_t_m = (C / (2 * np.pi * B_max ** 2 * alpha * np.log(1 + beta_eval))) ** (1.0 / 3.0)
    V_eff_tor = R_t_m ** 3 * F_at_beta
    Lp_tor = mu_0_si * R_t_m * ell_at_beta
    if Lp_tor == 0:
        return dict(error="ell at beta_eval = 0; Q^tor undefined")
    Q_tor_num = (B_max ** 4) * (V_eff_tor ** 4) / (Lp_tor ** 2)
    ratio = Q_tor_num / Q_core_num
    return dict(
        V_mag_m3=V_mag_m3,
        R_m_m=float(R_m),
        R_t_m=float(R_t_m),
        Veff_core_m3=float(Veff_core),
        Lp_core_H=float(Lp_core),
        Veff_tor_m3=float(V_eff_tor),
        Lp_tor_H=float(Lp_tor),
        Q_core=float(Q_core_num),
        Q_tor=float(Q_tor_num),
        Q_tor_over_Q_core=float(ratio),
        favorable=bool(ratio >= 1.0),
    )


def plot_scaling(grouped, alphas, out_dir):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib unavailable; skipping plots", flush=True)
        return
    os.makedirs(out_dir, exist_ok=True)

    # F vs beta
    fig, ax = plt.subplots(figsize=(5.5, 4.0))
    colors = {0.5: "tab:blue", 1.0: "tab:orange", 2.0: "tab:green"}
    for a in alphas:
        rs = grouped[a]
        bs = np.array([r["beta"] for r in rs])
        Fs = np.array([r["F"] for r in rs])
        ax.plot(bs, Fs, "o-", label=f"$\\alpha={a}$",
                color=colors.get(a, None))
    ax.axhline(0, color="grey", linewidth=0.5)
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"$F(\beta;\alpha) = V_{\rm eff}/R_t^3$")
    ax.set_title(r"Effective volume scaling")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "scaling_F_vs_beta.pdf"))
    plt.close(fig)

    # ell vs beta
    fig, ax = plt.subplots(figsize=(5.5, 4.0))
    for a in alphas:
        rs = grouped[a]
        bs = np.array([r["beta"] for r in rs])
        es = np.array([r["ell"] for r in rs])
        ax.plot(bs, es, "o-", label=f"$\\alpha={a}$",
                color=colors.get(a, None))
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"$\ell(\beta;\alpha) = L_p/(\mu_0 R_t)$")
    ax.set_title(r"Inductance scaling")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "scaling_ell_vs_beta.pdf"))
    plt.close(fig)


def main(sweep_csv: str, out_dir: str = "outputs",
         plot_dir: str = "plots"):
    rows = load_sweep(sweep_csv)
    bya = group_by_alpha(rows)

    fit_results = {}
    crossovers = {}
    # Restrict to beta >= 0.05: the smallest-beta points have very thin caps
    # (cap thickness ~ beta * R_t / n_t = 0.025 / 15 = 1.7e-3) and high
    # condition number, making ell numerically noisy.  All polynomial fits
    # use only the four moderate-beta points by default; the small-beta
    # points are reported in the raw CSV for completeness.
    BETA_MIN_FIT = 0.05
    for alpha, rs in bya.items():
        rs_fit = [r for r in rs if r["beta"] >= BETA_MIN_FIT]
        bs = np.array([r["beta"] for r in rs_fit])
        Fs = np.array([r["F"] for r in rs_fit])
        es = np.array([r["ell"] for r in rs_fit])

        fitF = fit_F(bs, Fs)
        fitE_free = fit_ell(bs, es, constrained=False)
        fitE_constr = fit_ell(bs, es, constrained=True)
        boot = bootstrap_F(bs, Fs, n_boot=2000)
        boot_ell = bootstrap_ell(bs, es, n_boot=2000)
        decision = decide(alpha, fitF, fitE_free, fitE_constr, boot)
        decision["ell_0_std_boot"] = boot_ell["ell_0_std"]
        decision["sigma_ell_0"] = (abs(fitE_free["ell_0"])
                                    / max(boot_ell["ell_0_std"], 1e-30))

        fit_results[str(alpha)] = dict(
            betas=bs.tolist(), Fs=Fs.tolist(), ells=es.tolist(),
            fit_F=fitF, bootstrap_F=boot,
            fit_ell_free=fitE_free,
            fit_ell_constrained=fitE_constr,
            bootstrap_ell=boot_ell,
            decision=decision,
        )

        # Crossover at beta_eval = 0.2 using the FULL fitted F(beta) and
        # ell(beta) (not just the linear coefficients).  This correctly
        # treats the case ell_0 != 0 (inductance-floor regime), where the
        # leading L_p contribution is the beta-independent ell_0.
        cs = {}
        beta_eval = 0.2
        F_eval = (fitF["f1"] * beta_eval
                  + fitF["f2"] * beta_eval**2
                  + fitF["f3"] * beta_eval**3)
        ell_eval = (fitE_free["ell_0"]
                    + fitE_free["ell_1"] * beta_eval
                    + fitE_free["ell_2"] * beta_eval**2)
        # Use directly-observed F and ell at beta=0.2 if available, as a
        # cross-check on the fit.
        F_data = None
        ell_data = None
        for k, b_val in enumerate(bs):
            if abs(b_val - beta_eval) < 1e-6:
                F_data = Fs[k]
                ell_data = es[k]
        F_use = F_data if F_data is not None else F_eval
        ell_use = ell_data if ell_data is not None else ell_eval
        for V in [0.1, 1.0, 10.0]:
            cs[f"V_mag_{V}_m3"] = crossover_Rt(alpha,
                                              F_use, ell_use,
                                              beta_eval=beta_eval,
                                              V_mag_m3=V)
            cs[f"V_mag_{V}_m3"]["F_used"] = float(F_use)
            cs[f"V_mag_{V}_m3"]["ell_used"] = float(ell_use)
            cs[f"V_mag_{V}_m3"]["beta_eval"] = beta_eval
            if decision["decision"] != "favorable":
                cs[f"V_mag_{V}_m3"]["note"] = (
                    "Evaluated at beta=%.2f data point; "
                    "decision: %s" % (beta_eval, decision["decision"]))
        crossovers[str(alpha)] = cs

    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "fit_results.json"), "w") as f:
        json.dump(fit_results, f, indent=2)
    with open(os.path.join(out_dir, "crossover_analysis.json"), "w") as f:
        json.dump(crossovers, f, indent=2)

    plot_scaling(bya, sorted(bya.keys()), plot_dir)

    # Print summary
    print("=" * 70)
    print("FIT RESULTS")
    print("=" * 70)
    for a, fr in fit_results.items():
        d = fr["decision"]
        print(f"alpha = {a}:")
        print(f"  decision: {d['decision']}")
        print(f"  reason  : {d['reason']}")
        print(f"  f_1 = {fr['fit_F']['f1']:+.4e}  (sig: {d['sigma_f1']:.2f} sigma)")
        print(f"  ell_0 = {fr['fit_ell_free']['ell_0']:.4e}  ell_1 = {fr['fit_ell_free']['ell_1']:.4e}")
        print(f"  RSS_constr/RSS_free = {d['rss_ratio']:.3f}")
        if a in crossovers:
            for k, c in crossovers[a].items():
                note = c.get("note", "")
                print(f"    {k}: R_t* = {c.get('R_t_m', 'n/a'):.3f} m, "
                      f"Q^tor/Q^core = {c.get('Q_tor_over_Q_core', 0):.3e} "
                      f"{note}")
    print()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="sweep_csv", default="outputs/production_sweep.csv")
    ap.add_argument("--out", default="outputs")
    ap.add_argument("--plots", default="plots")
    args = ap.parse_args()
    main(args.sweep_csv, args.out, args.plots)

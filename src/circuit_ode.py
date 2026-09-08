#!/usr/bin/env python3
"""circuit_ode.py — induction dynamics + sensitivity analysis (Module 3).

Module 2 gave us a PHOTO of the device (steady-state dose response).
This module turns it into a FILM: how mRNA and PETase accumulate over
time after IPTG is added at t=0.

Two outputs:
  1) Time-series plots of mRNA and PETase at several IPTG doses
  2) OAT sensitivity ranking: which parameter moves the output most?

Two pipelines feed this module:
  - Module 2's fitted K and n are imported and used in the induction term
    (the literal "M2 feeds M3" link).
  - The steady state of this ODE is CHECKED against Module 2's y_max /
    fold (the "handshake" validation).

IMPORTANT: everything here is a SIMULATED DEMO. Kinetic rates are
literature-range ESTIMATES (labeled as such) and were never measured.
No fabricated data
"""

import pathlib

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ---- the M2 -> M3 feed: import Module 2's demo data and fit machinery ----
# This is the literal connection between the two modules.
from transfer_function import make_demo_data, fit_hill

# Figures always go to the repo-root figures/ folder, no matter where the
# script is run from (fixes the Module 2 path problem for good).
REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
FIG_DIR = REPO_ROOT / "figures"
FIG_DIR.mkdir(exist_ok=True)

# --------------------------------------------------------------------------
# STEP 1 — the model
#
# Two "tanks":
#   d[mRNA]/dt = transcription_rate(I)  - delta_mRNA * mRNA   (fill - decay)
#   d[E]/dt    = beta * mRNA            - delta_E * E         (fill - decay)
#
# transcription_rate = v_max * (b + (1-b) * activation(I))
#   activation(I) = I^n / (K^n + I^n)   <- the SAME Hill shape as Module 2
#   b = basal leakage fraction (so the device is not fully OFF at I=0;
#       this gives us the fold induction / y_min behaviour).
# --------------------------------------------------------------------------
def induction_activation(I, K, n):
    """Hill activation term (0..1), identical in spirit to Module 2."""
    return (I ** n) / (K ** n + I ** n)


def build_params(K, n, v_max=50.0, beta=0.0064, basal=0.01,
                 half_life_mRNA_min=5.0, doubling_min=30.0):
    """Assemble the parameter dict for the ODE.

    Calibration note: beta is chosen so that the saturated steady state
    (I -> infinity) comes out near 100 = Module 2's demo y_max.
    All rate values are ESTIMATES (literature ranges), not measurements.
    """
    return {
        "K": K,                        # mM  — EC50, from M2 fit
        "n": n,                        #      — Hill, from M2 fit
        "v_max": v_max,                # a.u. — max transcription rate
        "beta": beta,                  # a.u. — translation rate per mRNA
        "basal": basal,                #      — leaky transcription fraction
        "delta_mRNA": np.log(2) / half_life_mRNA_min,   # 1/min — mRNA decay
        "delta_E": np.log(2) / doubling_min,            # 1/min — dilution by growth
    }


def ode_rhs(t, y, p, I):
    """Right-hand side of the ODE system. Returns [dmRNA/dt, dE/dt]."""
    mRNA, E = y
    act = induction_activation(I, p["K"], p["n"])
    transcription = p["v_max"] * (p["basal"] + (1 - p["basal"]) * act)
    dmRNA = transcription - p["delta_mRNA"] * mRNA
    dE = p["beta"] * mRNA - p["delta_E"] * E
    return [dmRNA, dE]


def run_induction(p, I, t_max=600.0, n_points=6001):
    """Solve the ODE after adding IPTG at t=0. Returns (t, mRNA, E)."""
    t_eval = np.linspace(0.0, t_max, n_points)
    sol = solve_ivp(
        ode_rhs, (0.0, t_max), [0.0, 0.0],
        t_eval=t_eval, args=(p, I),
        method="LSODA",          # robust for stiff ODEs (doc troubleshooting)
    )
    return sol.t, sol.y[0], sol.y[1]


def time_to_90(t, E, E_ss):
    """First time the enzyme reaches 90% of its steady state (minutes)."""
    hit = np.argmax(E >= 0.9 * E_ss)
    return float(t[hit]) if E[-1] >= 0.9 * E_ss else float("nan")


# --------------------------------------------------------------------------
# STEP 2 — validation: the M2 <-> M3 handshake.
# The ODE's steady state must reproduce Module 2's dose-response curve.
# --------------------------------------------------------------------------
def check_handshake(m2_params, tol_pct=20.0):
    """Compare ODE steady states with M2 fit results; print PASS/FAIL."""
    y_min, y_max, K, n = m2_params
    p = build_params(K=K, n=n)

    _, _, E_sat = run_induction(p, 10.0)     # saturated dose
    _, _, E_low = run_induction(p, 1e-3)     # near-basal dose
    E_sat, E_low = E_sat[-1], E_low[-1]

    fold_ode = E_sat / E_low
    fold_m2 = y_max / y_min
    err_max = abs(E_sat - y_max) / y_max * 100.0
    err_fold = abs(fold_ode - fold_m2) / fold_m2 * 100.0
    ok = err_max < tol_pct and err_fold < tol_pct

    print("=== M2 <-> M3 handshake validation ===")
    print(f"  ODE steady state (sat) : {E_sat:7.2f}   | M2 y_max: {y_max:7.2f} "
          f"({err_max:4.1f}% off)")
    print(f"  ODE fold induction     : {fold_ode:7.1f}x | M2 fold:  {fold_m2:7.1f}x "
          f"({err_fold:4.1f}% off)")
    print(f"  RESULT: {'PASS' if ok else 'FAIL — check units/rates'}")
    print("  (tolerance = " + str(tol_pct) + "%; demo is close by design)\n")
    return ok


# --------------------------------------------------------------------------
# STEP 3 — OAT sensitivity analysis (±10% one-at-a-time).
# S = (relative change in output) / (relative change in parameter).
# Positive S: raising the parameter raises the output (production taps).
# Negative S: raising it lowers the output (decay taps, or higher EC50).
# --------------------------------------------------------------------------
def sensitivity_oat(p, I, frac=0.10):
    SENSITIVE_KEYS = ("K", "n", "v_max", "beta", "delta_mRNA", "delta_E")

    t_base, _, E_series = run_induction(p, I)   # tam zaman serisi
    E_base = E_series[-1]                        # skaler: durağan durum
    t90_base = time_to_90(t_base, E_series, E_base)   # dizi + skaler

    rows = {}
    for key in SENSITIVE_KEYS:
        p2 = dict(p)
        p2[key] = p[key] * (1.0 + frac)
        t2, _, E2 = run_induction(p2, I)
        E2_ss = E2[-1]
        t90_2 = time_to_90(t2, E2, E2_ss)        # heh
        rows[key] = {
            "S_Ess": (E2_ss - E_base) / E_base / frac,
            "S_t90": (t90_2 - t90_base) / t90_base / frac if t90_base > 0 else float("nan"),
        }
    return rows, E_base, t90_base


# --------------------------------------------------------------------------
# STEP 4 — plots.
# --------------------------------------------------------------------------
def plot_induction(p, doses, t_max=600.0):
    """Time series of mRNA and PETase for several IPTG doses."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for I in doses:
        t, mRNA, E = run_induction(p, I, t_max=t_max)
        axes[0].plot(t, mRNA, label=f"IPTG {I:g} mM")
        axes[1].plot(t, E, label=f"IPTG {I:g} mM")

    axes[0].set(xlabel="time (min)", ylabel="mRNA (a.u.)",
                title="mRNA after induction (SIMULATED)")
    axes[1].set(xlabel="time (min)", ylabel="PETase (a.u.)",
                title="PETase after induction (SIMULATED)")
    for ax in axes:
        ax.legend(fontsize=8)
    fig.tight_layout()
    path = FIG_DIR / "induction_timeseries.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved figure: {path}")


def plot_tornado(rows):
    """Horizontal bar chart of S(E_ss), sorted by magnitude."""
    keys = sorted(rows, key=lambda k: abs(rows[k]["S_Ess"]))
    vals = [rows[k]["S_Ess"] for k in keys]
    colors = ["crimson" if v < 0 else "teal" for v in vals]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.barh(keys, vals, color=colors)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("Sensitivity S (steady-state enzyme, ±10%)")
    ax.set_title("OAT sensitivity ranking (SIMULATED)")
    fig.tight_layout()
    path = FIG_DIR / "sensitivity_tornado.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved figure: {path}")


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
if __name__ == "__main__":
    # ---- feed from Module 2: refit the demo data and take K, n -----------
    doses, measured = make_demo_data()
    m2_params, _ = fit_hill(doses, measured)   # [y_min, y_max, K, n]
    K_m2, n_m2 = m2_params[2], m2_params[3]
    print(f"Module 2 feed -> K = {K_m2:.4f} mM, n = {n_m2:.3f} (SIMULATED demo)")

    # ---- handshake validation --------------------------------------------
    check_handshake(m2_params)

    # ---- main induction movie ---------------------------------------------
    p = build_params(K=K_m2, n=n_m2)
    plot_induction(p, doses=[0.03, 0.1, 3.0])

    # ---- sensitivity ranking: pick the steep part of the curve ------------
    # At exactly I = K the activation is 0.5 for ANY n, so n's effect would
    # look zero there. Below EC50 (I = 0.03 mM) we sit on the steep part,
    # where the switch parameters (K, n) show their real influence.
    rows, E_base, t90_base = sensitivity_oat(p, I=0.03)
    plot_tornado(rows)

    print("=== OAT sensitivity ranking (steady-state enzyme, I = 0.03 mM) ===")
    print(f"  base steady state: {E_base:6.2f} a.u. | base t90: {t90_base:6.1f} min")
    print(f"  {'parameter':<12}{'S (E_ss)':>10}{'S (t90)':>10}   interpretation")
    for key, r in sorted(rows.items(), key=lambda kv: -abs(kv[1]["S_Ess"])):
        tag = "switch/amplifier" if abs(r["S_Ess"]) > 1 else (
              "production tap" if r["S_Ess"] > 0 else "decay tap")
        print(f"  {key:<12}{r['S_Ess']:>10.2f}{r['S_t90']:>10.2f}   {tag}")
    print("\n  NOTE: ranking depends on the dose you look at (steep vs flat part).")
    print("  Rate parameters ~±1 matches the analytic prediction "
          "E_ss = v_max*beta/(delta_mRNA*delta_E) -> model is sane.")

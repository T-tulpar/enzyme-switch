#!/usr/bin/env python3
"""transfer_function.py

Answers the question: "how much enzyme do I get for each dose of IPTG?"
We fit a Hill curve to dose-response data and read off 3 numbers:
  EC50 (K), Hill coefficient (n), fold induction (dynamic range).

NOTE: the data used here is SIMULATED DEMO DATA (no wet-lab measurement).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# ----------------------------------------------------------------------
# THE HILL EQUATION — the S-curve we expect from this circuit.
#   y(I) = y_min + (y_max - y_min) * (I**n) / (K**n + I**n)
#
#   y     = enzyme output
#   I     = IPTG dose (input)
#   K     = EC50, the dose giving half-max output (switch sensitivity)
#   n     = Hill coefficient, curve steepness (n>1 = switchy/digital)
#   y_min = basal (leaky) output at zero inducer
#   y_max = saturated output at high inducer
# ----------------------------------------------------------------------
def hill_activation(I, y_min, y_max, K, n):
    """Evaluate the Hill activation curve at dose I (vectorised)."""
    return y_min + (y_max - y_min) * (I**n) / (K**n + I**n)

# ----------------------------------------------------------------------
# STEP 1 — build SIMULATED demo data (not a measurement!).
#
# True underlying parameters (the "answer key" of the simulation):
#   K ~ 0.1 mM IPTG, n ~ 1.8, fold induction ~100x.
# These are biologically plausible for LacI/IPTG/T7, but they are DEMO.
# We generate fake noisy measurements around this true curve, so that
# the fit later should recover values close to these truth parameters.
# ----------------------------------------------------------------------
def make_demo_data(rng_seed=42):
    """Return (doses, measurements) simulated around known Hill params."""
    # hidden truth parameters — only used to CREATE the fake data
    true_K = 0.1    # mM IPTG (EC50)
    true_n = 1.8    # Hill coefficient
    true_ymin = 1.0 # basal output (relative units)
    true_ymax = 100.0 # saturated output (100x above basal)

    # doses spread across a log range (0.001 to 10 mM IPTG), 40 points
    doses = np.logspace(-3, 1, 40)

    # clean signal from the true curve
    signal = hill_activation(doses, true_ymin, true_ymax, true_K, true_n)

    # add random measurement noise (proportional noise, biological style)
    rng = np.random.default_rng(rng_seed)
    noise = rng.normal(0, 0.05, size=doses.size)  # ~5% noise
    measured = signal * (1 + noise)

    # clip at zero (can't have negative enzyme)
    measured = np.maximum(measured, 0)

    return doses, measured

# ----------------------------------------------------------------------
# STEP 2 — fit the Hill curve to the data.
# curve_fit tweaks (y_min, y_max, K, n) until its curve best matches data.
# We give log-space-friendly initial guesses + bounds to help it converge.
# ----------------------------------------------------------------------
def fit_hill(doses, measured):
    """Fit the Hill activation model. Return fitted params + covariance."""
    # initial guesses: K=0.5, n=1, y_min=min(data), y_max=max(data)
    p0 = [measured.min(), measured.max(), 0.5, 1.0]
    # bounds keep parameters physically sensible
    bounds = (
        [0,   0,   1e-4, 0.1],   # lower bounds
        [1e3, 1e4, 10,   8],     # upper bounds
    )
    params, cov = curve_fit(
        hill_activation, doses, measured, p0=p0, bounds=bounds
    )
    return params, cov

# ----------------------------------------------------------------------
# STEP 3 — validation: are the residuals small? Is the fit good?
# We also package the "device ID card" (the 3 headline numbers).
# ----------------------------------------------------------------------
def summarize(params, doses, measured):
    """Print the device summary and return a dict of headline numbers."""
    y_min, y_max, K, n = params

    # fit quality: compare model vs data
    predicted = hill_activation(doses, *params)
    residuals = measured - predicted
    rms_resid = float(np.sqrt(np.mean(residuals**2)))

    fold = y_max / y_min  # dynamic range = fold induction

    print("=== Transfer function fit (SIMULATED DEMO DATA) ===")
    print(f"  EC50 (K)          : {K:8.4f} mM IPTG")
    print(f"  Hill coefficient n: {n:8.3f}")
    print(f"  Basal (y_min)     : {y_min:8.2f}")
    print(f"  Saturated (y_max) : {y_max:8.2f}")
    print(f"  Fold induction    : {fold:8.1f}x")
    print(f"  Fit RMS residual  : {rms_resid:8.4f} (small = good fit)")
    print("\n  NOTE: this is SIMULATED demo data, not a wet-lab measurement.")

    return {"EC50_mM": K, "n": n, "ymin": y_min, "ymax": y_max, "fold": fold}

# ----------------------------------------------------------------------
# STEP 4 — plot data points + fitted curve -> figures/transfer_function.png
# ----------------------------------------------------------------------
def plot_fit(doses, measured, params):
    """Save the dose-response figure to figures/transfer_function.png."""
    # dense curve for a smooth line
    I_smooth = np.logspace(-3, 1, 200)
    y_smooth = hill_activation(I_smooth, *params)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(doses, measured, label="data (SIMULATED)", color="gray", s=30)
    ax.plot(I_smooth, y_smooth, label="Hill fit", color="crimson")
    ax.set_xscale("log")            # IPTG dose is better shown on log axis
    ax.set_xlabel("IPTG dose (mM)")
    ax.set_ylabel("Enzyme output (relative)")
    ax.set_title("EnzymeSwitch — PETase transfer function (SIMULATED)")
    ax.legend()
    fig.tight_layout()
    fig.savefig("figures/transfer_function.png", dpi=150)
    plt.close(fig)                  # don't leave a hanging window
    print("  Saved figure: figures/transfer_function.png")

# ----------------------------------------------------------------------
# main — only runs when you type: python src/transfer_function.py
# ----------------------------------------------------------------------
if __name__ == "__main__":
    doses, measured = make_demo_data()
    params, cov = fit_hill(doses, measured)
    summary = summarize(params, doses, measured)
    plot_fit(doses, measured, params)

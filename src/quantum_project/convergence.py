"""Numerical error metrics for probability trajectories on matching time grids."""

import numpy as np

from .states import PSI_PLUS


def sweep_probabilities(states):
    """Return columns P00, P0r, Pr0, Prr, P_Psi+ without clipping."""
    states = np.asarray(states, dtype=complex)
    if states.ndim != 2 or states.shape[1] != 4 or not np.all(np.isfinite(states)):
        raise ValueError("states must be a finite array of shape (n, 4)")
    return np.column_stack((np.abs(states)**2, np.abs(states @ PSI_PLUS.conj())**2))


def significant_tolerance(value, digits=4):
    """Half a unit in the last significant place at a positive reference scale.

    This is an absolute tolerance, not a fixed number of decimal places.
    A zero reference has no significant scale and is rejected.
    """
    if not isinstance(digits, int) or isinstance(digits, bool) or digits < 1:
        raise ValueError("digits must be a positive integer")
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise ValueError("reference scale must be finite and positive")
    return float(0.5 * 10.0**(np.floor(np.log10(value)) - digits + 1))


def compare_probabilities(actual, reference, reference_peak_fidelity=None):
    """Compare (n, 5) probabilities at identical times using four-place scales.

    Trajectory tolerances use each reference observable's peak scale, avoiding
    undefined relative errors at zeros. Final and peak Bell fidelities also
    require equal four-significant-digit rounding. The time of the sampled
    maximum is reported elsewhere and is not certified by this criterion.
    """
    actual, reference = np.asarray(actual, dtype=float), np.asarray(reference, dtype=float)
    if (actual.shape != reference.shape or actual.ndim != 2 or actual.shape[1] != 5
            or len(actual) == 0 or not np.all(np.isfinite(actual))
            or not np.all(np.isfinite(reference))):
        raise ValueError("probabilities must be finite, matching nonempty (n, 5) arrays")
    error_00 = float(np.max(np.abs(actual[:, 0] - reference[:, 0])))
    error_bell = float(np.max(np.abs(actual[:, 4] - reference[:, 4])))
    final = float(actual[-1, 4])
    final_ref = float(reference[-1, 4])
    peak = float(actual[:, 4].max())
    peak_ref = float(reference[:, 4].max())
    if reference_peak_fidelity is not None:
        supplied_peak = float(reference_peak_fidelity)
        if not np.isfinite(supplied_peak) or supplied_peak < peak_ref - 1e-12:
            raise ValueError("reference peak must be finite and at least the sampled peak")
        peak_ref = supplied_peak
    tolerances = {
        "P_00": significant_tolerance(float(reference[:, 0].max())),
        "P_Psi_plus": significant_tolerance(peak_ref),
        "final_fidelity": significant_tolerance(final_ref),
        "peak_fidelity": significant_tolerance(peak_ref),
    }
    final_error, peak_error = abs(final - final_ref), abs(peak - peak_ref)
    ratios = [error_00 / tolerances["P_00"], error_bell / tolerances["P_Psi_plus"],
              final_error / tolerances["final_fidelity"], peak_error / tolerances["peak_fidelity"]]
    rounded_equal = format(final, ".4g") == format(final_ref, ".4g") and format(peak, ".4g") == format(peak_ref, ".4g")
    return {
        "final_bell_fidelity": final, "reference_final_bell_fidelity": final_ref,
        "final_fidelity_error": final_error, "maximum_P_00_error": error_00,
        "maximum_P_Psi_plus_error": error_bell,
        "maximum_bell_fidelity_difference": error_bell,
        "peak_bell_fidelity": peak, "peak_fidelity_error": peak_error,
        "tolerances": tolerances, "maximum_error_ratio": float(max(ratios)),
        "scalar_rounding_matches": rounded_equal,
        "passes": bool(max(ratios) <= 1 and rounded_equal),
    }


def largest_passing_timestep(rows):
    """Largest tested dt whose row and all finer tested rows pass.

    Return None if no such row exists. This avoids accepting an isolated
    cancellation at a coarse step; no claim about untested steps is made.
    """
    if not rows:
        raise ValueError("at least one convergence row is required")
    ordered = sorted(rows, key=lambda row: row["dt_s"])
    dts = np.array([row["dt_s"] for row in ordered])
    if not np.all(np.isfinite(dts)) or np.any(dts <= 0) or len(np.unique(dts)) != len(dts):
        raise ValueError("timesteps must be unique, finite and positive")
    best = None
    for row in ordered:
        if not row["passes"]:
            break
        best = row["dt_s"]
    return best

"""Generate the static Rydberg figure, raw data and analytical residual report.

Run from the repository root after installing the project:
    python scripts/figure_5_1.py
"""

import json
from pathlib import Path
import platform

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from quantum_project import KET_00, evolve_static, rydberg_hamiltonian
from quantum_project.convergence import sweep_probabilities as probabilities
from quantum_project.static_limits import blockade_state, noninteracting_state
from quantum_project.units import mhz_to_angular_frequency, microseconds_to_seconds


ROOT = Path(__file__).resolve().parents[1]
OMEGA = mhz_to_angular_frequency(1.0)
PLOT_V_MHZ = (0.0, 0.5, 1.0, 5.0, 100.0)
VALIDATION_CASES = (
    ("noninteracting", 0.0), ("noninteracting", 0.001), ("noninteracting", 0.01),
    ("blockade", 10.0), ("blockade", 30.0), ("blockade", 100.0),
    ("blockade", -100.0),
)
OBSERVABLES = ("P_00", "P_0r", "P_r0", "P_rr", "P_Psi_plus")
TOLERANCE = 1e-12


def physical_checks(values):
    """Measure and enforce probability conservation, exchange symmetry and bounds."""
    if not np.all(np.isfinite(values)):
        raise ValueError("nonfinite probabilities")
    checks = {
        "max_probability_normalisation_error": float(np.max(np.abs(values[:, :4].sum(axis=1) - 1))),
        "max_exchange_symmetry_error": float(np.max(np.abs(values[:, 1] - values[:, 2]))),
        "minimum_probability": float(values.min()),
        "maximum_probability": float(values.max()),
    }
    if (checks["max_probability_normalisation_error"] > TOLERANCE
            or checks["max_exchange_symmetry_error"] > TOLERANCE
            or checks["minimum_probability"] < -TOLERANCE
            or checks["maximum_probability"] > 1 + TOLERANCE):
        raise ValueError(f"physical check failed: {checks}")
    return checks


def simulate(times, interaction_mhz):
    hbar = rydberg_hamiltonian(OMEGA, 0.0, 0.0, mhz_to_angular_frequency(interaction_mhz))
    states = evolve_static(KET_00, times, hbar)
    values = probabilities(states)
    return states, values, physical_checks(values)


def residual_metrics(residual):
    return {
        name: {"max_absolute": float(np.max(np.abs(residual[:, index]))),
               "rms": float(np.sqrt(np.mean(residual[:, index]**2)))}
        for index, name in enumerate(OBSERVABLES)
    }


def save_figure(figure, path):
    figure.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(figure)


def plot_dynamics(times_us, values, directory):
    colors = ("#0072B2", "#ef8a62", "#d6604d", "#b2182b", "#67001f")
    styles = ("-", "--", "-.", ":", (0, (5, 1.5)))
    fig, axes = plt.subplots(2, 1, figsize=(7.1, 6.4), sharex=True, layout="constrained")
    for interaction, case, color, style in zip(PLOT_V_MHZ, values, colors, styles):
        for ax, index in zip(axes, (0, 4)):
            ax.plot(times_us, case[:, index], color=color, linestyle=style,
                    label=f"{interaction:g}", linewidth=1.7)
    for ax, label in zip(axes, (r"(a) $P_{00}$", r"(b) $P_{\Psi^+}$")):
        ax.set(ylabel="Probability", ylim=(-0.025, 1.025), xlim=(0, 1))
        ax.set_yticks(np.linspace(0, 1, 5))
        ax.text(0.025, 0.88, label, transform=ax.transAxes,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85})
        ax.grid(alpha=0.2, linewidth=0.6)
    axes[0].legend(title=r"$V/(2\pi)$ (MHz)", ncol=5, loc="lower left",
                   bbox_to_anchor=(0, 1.01, 1, 0.15), mode="expand", borderaxespad=0)
    axes[-1].set_xlabel("Time (microseconds)")
    fig.suptitle(r"Two-atom Rydberg dynamics: $\Omega/(2\pi)=1$ MHz, $\Delta=\phi=0$",
                 fontsize=12)
    save_figure(fig, directory / "figure_5_1")


def plot_residuals(times_us, residuals, directory):
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.0), sharex=True, layout="constrained")
    # Exact-zero and finite-blockade residuals have very different scales.
    for column, case_index in enumerate((0, 5)):
        for row, observable in enumerate((0, 4)):
            ax = axes[row, column]
            ax.plot(times_us, residuals[case_index, :, observable], color="#0072B2", linewidth=1)
            ax.axhline(0, color="0.4", linewidth=0.6)
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
            ax.set_xlim(0, 1)
            ax.grid(alpha=0.2, linewidth=0.6)
            ax.set_ylabel((r"$P_{00}$" if row == 0 else r"$P_{\Psi^+}$") + " residual")
    axes[0, 0].set_title(r"$V=0$: independent atoms")
    axes[0, 1].set_title(r"$V/(2\pi)=100$ MHz: ideal blockade")
    for ax in axes[-1]:
        ax.set_xlabel("Time (microseconds)")
    fig.suptitle("Probability residuals: numerical − analytical", fontsize=12)
    save_figure(fig, directory / "static_residuals")


def main():
    figures, data = ROOT / "results/figures", ROOT / "results/data"
    figures.mkdir(parents=True, exist_ok=True)
    data.mkdir(parents=True, exist_ok=True)
    times = np.linspace(0, microseconds_to_seconds(1), 4001)
    times_us = times / microseconds_to_seconds(1)
    report = {
        "source_interpretation": "Source specifies Omega/(2pi)=1 MHz, blue V=0, red curves spanning 0–100 MHz; individual red-curve V values are not enumerated. The five plot values are repository numerical sampling choices, not an exact source curve list or colour density.",
        "parameters": {"omega_rad_s": OMEGA, "delta_rad_s": 0.0, "phase_rad": 0.0,
                       "initial_state": "|00>", "basis_order": ["00", "0r", "r0", "rr"],
                       "plot_V_over_2pi_MHz": PLOT_V_MHZ, "time_start_s": float(times[0]),
                       "time_end_s": float(times[-1]), "time_samples": len(times),
                       "physical_check_tolerance": TOLERANCE},
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "matplotlib": matplotlib.__version__},
        "residual_definition": "numerical probability minus limiting analytical probability",
        "rms_definition": "sqrt(mean(residual**2)) over the uniformly sampled inclusive time grid",
        "plot_checks": [], "validation": [],
    }
    plot_states, plot_values = [], []
    for interaction in PLOT_V_MHZ:
        states, values, checks = simulate(times, interaction)
        plot_states.append(states)
        plot_values.append(values)
        report["plot_checks"].append({"V_over_2pi_MHz": interaction, **checks})

    validation_states, analytical_states, residuals = [], [], []
    fine_times = np.linspace(times[0], times[-1], 8001)
    for limit, interaction in VALIDATION_CASES:
        analytical = noninteracting_state if limit == "noninteracting" else blockade_state
        states, values, checks = simulate(times, interaction)
        expected = analytical(times, OMEGA)
        residual = values - probabilities(expected)
        metrics = residual_metrics(residual)
        # Double the sampling to check that narrow finite-V oscillations are resolved.
        _, fine_values, fine_checks = simulate(fine_times, interaction)
        fine_metrics = residual_metrics(fine_values - probabilities(analytical(fine_times, OMEGA)))
        report["validation"].append({
            "limit": limit, "V_over_2pi_MHz": interaction,
            "residuals": metrics, "physical_checks": checks,
            "refined_grid_samples": len(fine_times), "refined_residuals": fine_metrics,
            "refined_physical_checks": fine_checks,
        })
        validation_states.append(states)
        analytical_states.append(expected)
        residuals.append(residual)

    plot_values, residuals = np.array(plot_values), np.array(residuals)
    all_checks = report["plot_checks"] + [
        case[key] for case in report["validation"]
        for key in ("physical_checks", "refined_physical_checks")
    ]
    report["overall_checks"] = {
        "max_probability_normalisation_error": max(c["max_probability_normalisation_error"] for c in all_checks),
        "max_exchange_symmetry_error": max(c["max_exchange_symmetry_error"] for c in all_checks),
        "minimum_probability": min(c["minimum_probability"] for c in all_checks),
        "maximum_probability": max(c["maximum_probability"] for c in all_checks),
    }
    np.savez_compressed(
        data / "figure_5_1.npz", times_s=times, plot_V_over_2pi_MHz=PLOT_V_MHZ,
        omega_rad_s=OMEGA, delta_rad_s=0.0, phase_rad=0.0,
        basis_labels=np.array(["00", "0r", "r0", "rr"]),
        observable_labels=np.array(OBSERVABLES), plot_states=np.array(plot_states),
        plot_probabilities=plot_values,
        validation_V_over_2pi_MHz=[v for _, v in VALIDATION_CASES],
        validation_limits=np.array([limit for limit, _ in VALIDATION_CASES]),
        validation_states=np.array(validation_states), analytical_states=np.array(analytical_states),
        validation_probabilities=np.array([probabilities(s) for s in validation_states]),
        analytical_probabilities=np.array([probabilities(s) for s in analytical_states]),
        residuals=residuals,
    )
    (data / "static_validation.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42}):
        plot_dynamics(times_us, plot_values, figures)
        plot_residuals(times_us, residuals, figures)
    print("Limit          V/(2pi) MHz  Observable    Max absolute       RMS")
    for case in report["validation"]:
        for observable in ("P_00", "P_Psi_plus"):
            metric = case["residuals"][observable]
            print(f"{case['limit']:14s} {case['V_over_2pi_MHz']:11g}  {observable:11s} "
                  f"{metric['max_absolute']:.8e}  {metric['rms']:.8e}")
    print(json.dumps(report["overall_checks"], indent=2))
    print(f"Figures: {figures}\nData and report: {data}")


if __name__ == "__main__":
    main()

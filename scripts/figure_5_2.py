"""Generate detuning-sweep figures and convergence results for explicit parameters.

Milestone reproduction (documented sampling of the source's 0–100 MHz range):
    python scripts/figure_5_2.py --milestone --interaction-mhz 0 0.5 1 5 100 --phase 0
"""

import argparse
import csv
import json
from pathlib import Path
import platform

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

from quantum_project import KET_00, evolve_trotter, rydberg_hamiltonian
from quantum_project.convergence import compare_probabilities, largest_passing_timestep, sweep_probabilities
from quantum_project.sweeps import cubic_detuning, linear_detuning
from quantum_project.units import mhz_to_angular_frequency


ROOT = Path(__file__).resolve().parents[1]
TMAX = 1e-6
OMEGA = mhz_to_angular_frequency(1)
DELTA_MAX = mhz_to_angular_frequency(50)
PROFILES = {"cubic": cubic_detuning, "linear": linear_detuning}
STEP_COUNTS = (131072, 65536, 32768, 16384, 8192, 4096, 2048, 1024, 512, 256)
# Fine-reference and audit settings are documented in docs/sweep_validation.md.
AUDIT_SAMPLES = 40001
REFERENCE_BASE_STEPS = 2000
REFERENCE_PROBABILITY_TOLERANCE = 1e-9
PHYSICAL_TOLERANCE = 1e-10  # Allows roundoff over >100,000 unitary steps.


def physical_checks(states):
    values = sweep_probabilities(states)
    checks = {
        "maximum_norm_error": float(np.max(np.abs(np.linalg.norm(states, axis=1) - 1))),
        "maximum_probability_sum_error": float(np.max(np.abs(values[:, :4].sum(axis=1) - 1))),
        "maximum_exchange_symmetry_error": float(np.max(np.abs(values[:, 1] - values[:, 2]))),
        "minimum_probability": float(values.min()), "maximum_probability": float(values.max()),
    }
    if (checks["maximum_norm_error"] > PHYSICAL_TOLERANCE
            or checks["maximum_probability_sum_error"] > PHYSICAL_TOLERANCE
            or checks["maximum_exchange_symmetry_error"] > PHYSICAL_TOLERANCE
            or checks["minimum_probability"] < -PHYSICAL_TOLERANCE
            or checks["maximum_probability"] > 1 + PHYSICAL_TOLERANCE):
        raise ValueError(f"Physical check failed: {checks}")
    return checks


def reference_solution(callback, refinement):
    solution = solve_ivp(
        lambda t, psi: -1j * callback(t) @ psi, (0, TMAX), KET_00.copy(),
        method="DOP853", dense_output=True, rtol=1e-11 if refinement == 1 else 1e-13,
        atol=1e-13 if refinement == 1 else 1e-15, max_step=TMAX / (REFERENCE_BASE_STEPS * refinement),
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution


def reference_peak(reference):
    times = np.linspace(0, TMAX, AUDIT_SAMPLES)
    values = sweep_probabilities(reference.sol(times).T)[:, 4]
    index = int(np.argmax(values))
    candidates = [(float(values[0]), 0.0), (float(values[-1]), TMAX), (float(values[index]), float(times[index]))]
    if 0 < index < len(times) - 1:
        # Optimise in microseconds to keep the bounded minimiser well scaled.
        fit = minimize_scalar(
            lambda time_us: -float(sweep_probabilities(reference.sol([time_us * 1e-6]).T)[0, 4]),
            bounds=(times[index - 1] * 1e6, times[index + 1] * 1e6),
            method="bounded", options={"xatol": 1e-12},
        )
        if fit.success:
            candidates.append((-float(fit.fun), float(fit.x * 1e-6)))
    return max(candidates)


def run_case(profile, interaction_mhz, phase, step_counts, family_check=False):
    sweep = PROFILES[profile]
    callback = lambda t: rydberg_hamiltonian(OMEGA, sweep(t, DELTA_MAX, TMAX), phase,
                                            mhz_to_angular_frequency(interaction_mhz))
    reference = reference_solution(callback, 1)
    refined = reference_solution(callback, 2)
    audit_times = np.linspace(0, TMAX, AUDIT_SAMPLES)
    audit_states = refined.sol(audit_times).T
    reference_difference = float(np.max(np.abs(
        sweep_probabilities(reference.sol(audit_times).T) - sweep_probabilities(audit_states))))
    if reference_difference > REFERENCE_PROBABILITY_TOLERANCE:
        raise ValueError("Reference refinement did not reach 1e-9 probability agreement")
    peak, peak_time = reference_peak(refined)
    final = float(sweep_probabilities(refined.sol([TMAX]).T)[0, 4])
    rows, saved = [], {}
    finer_values = None
    for steps in sorted(step_counts, reverse=True):
        times = np.linspace(0, TMAX, steps + 1)
        states = evolve_trotter(KET_00, times, callback)
        expected = refined.sol(times).T
        values = sweep_probabilities(states)
        metrics = compare_probabilities(values, sweep_probabilities(expected), reference_peak_fidelity=peak)
        if family_check:
            # Plot stability is separate from the benchmark's scalar rounding criterion.
            metrics["benchmark_criterion_passes"] = metrics["passes"]
            metrics["passes"] = max(metrics["maximum_P_00_error"],
                                   metrics["maximum_P_Psi_plus_error"]) <= 5e-5
        row = {"profile": profile, "V_over_2pi_MHz": interaction_mhz, "steps": steps,
               "dt_s": TMAX / steps, **metrics, "physical_checks": physical_checks(states),
               "validation_scope": "family plot stability" if family_check else "full convergence benchmark"}
        index = int(np.argmax(values[:, 4]))
        row["sampled_peak_time_s"] = float(times[index])
        row["reference_peak_time_s"] = peak_time
        if family_check and finer_values is not None:
            difference = float(np.max(np.abs(values[:, [0, 4]] - finer_values[::2, [0, 4]])))
            row["halved_dt_maximum_plot_difference"] = difference
            row["passes"] = row["passes"] and difference <= 5e-5
        finer_values = values
        rows.append(row)
        # Retain the finest result, then replace it only by a larger passing step.
        if not saved or largest_passing_timestep(rows) == row["dt_s"]:
            saved = {"times_s": times, "states": states, "probabilities": values,
                     "reference_states": expected, "reference_probabilities": sweep_probabilities(expected),
                     "detuning_rad_s": sweep(times, DELTA_MAX, TMAX)}
        print(f"{profile:6s} V/(2pi)={interaction_mhz:g} MHz dt={row['dt_s']*1e9:.6g} ns "
              f"error/tolerance={row['maximum_error_ratio']:.4g} pass={row['passes']}", flush=True)
    largest = largest_passing_timestep(rows)
    summary = {
        "profile": profile, "V_over_2pi_MHz": interaction_mhz,
        "reference_final_fidelity": final, "reference_maximum_fidelity": peak,
        "reference_time_of_maximum_s": peak_time,
        "reference_refinement_max_probability_difference": reference_difference,
        "reference_physical_checks": physical_checks(audit_states),
        "largest_tested_passing_dt_s": largest,
        "selected_steps": len(saved["times_s"]) - 1,
        "selected_final_fidelity": float(saved["probabilities"][-1, 4]),
        "selected_maximum_fidelity": float(saved["probabilities"][:, 4].max()),
        "selected_sampled_peak_time_s": float(saved["times_s"][np.argmax(saved["probabilities"][:, 4])]),
        "validation_scope": "family plot stability" if family_check else "full convergence benchmark",
    }
    return rows, summary, saved


def save_figure(fig, path):
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_figures(cases, summaries, rows, figures, case_label="User-specified parameters"):
    colors = ("#0072B2", "#ef8a62", "#d6604d", "#b2182b", "#67001f")
    interactions = list(dict.fromkeys(summary["V_over_2pi_MHz"] for summary in summaries))
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 8), sharex=True, layout="constrained")
    times = np.linspace(0, TMAX, 2001)
    axes[0].plot(times * 1e6, cubic_detuning(times, DELTA_MAX, TMAX) / (2 * np.pi * 1e6), color="#333333")
    axes[0].set_ylabel(r"$\Delta/(2\pi)$ (MHz)")
    for index, interaction in enumerate(interactions):
        case = cases[("cubic", interaction)]
        for ax, column in zip(axes[1:], (0, 4)):
            ax.plot(case["times_s"] * 1e6, case["probabilities"][:, column],
                    color=colors[index % len(colors)], label=f"{interaction:g} MHz", lw=1.5)
    for ax, label in zip(axes[1:], (r"$P_{00}$", r"$P_{\Psi^+}$")):
        ax.set_ylabel(label)
        ax.set_ylim(-0.025, 1.025)
    axes[1].legend(title=r"$V/(2\pi)$", ncol=min(3, len(interactions)), fontsize=9)
    for ax in axes:
        ax.set_xlim(0, 1)
        ax.grid(alpha=0.2)
    axes[-1].set_xlabel("Time (microseconds)")
    fig.suptitle(case_label + "\n" + r"Cubic sweep: $\Omega/(2\pi)=1$ MHz, $\Delta_{\max}/(2\pi)=50$ MHz")
    save_figure(fig, figures / "figure_5_2")

    fig, axes = plt.subplots(3, 1, figsize=(7.2, 8), sharex=True, layout="constrained")
    for profile, style in (("cubic", "-"), ("linear", "--")):
        axes[0].plot(times * 1e6, PROFILES[profile](times, DELTA_MAX, TMAX) / (2 * np.pi * 1e6),
                     linestyle=style, label=profile)
        for index, interaction in enumerate(interactions):
            case = cases[(profile, interaction)]
            for ax, column in zip(axes[1:], (0, 4)):
                ax.plot(case["times_s"] * 1e6, case["probabilities"][:, column],
                        linestyle=style, color=colors[index % len(colors)], lw=1.4,
                        label=f"{profile}, {interaction:g} MHz")
    axes[0].set_ylabel(r"$\Delta/(2\pi)$ (MHz)")
    axes[0].legend()
    for ax, label in zip(axes[1:], (r"$P_{00}$", r"$P_{\Psi^+}$")):
        ax.set_ylabel(label)
        ax.set_ylim(-0.025, 1.025)
    axes[1].legend(fontsize=8, ncol=2)
    for ax in axes:
        ax.grid(alpha=0.2)
        ax.set_xlim(0, 1)
    axes[-1].set_xlabel("Time (microseconds)")
    fig.suptitle(case_label + "\nCubic and linear detuning: identical endpoints and duration")
    save_figure(fig, figures / "sweep_comparison")

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.2), sharex=True, layout="constrained")
    for profile, style in (("cubic", "-"), ("linear", "--")):
        for index, interaction in enumerate(interactions):
            subset = sorted((row for row in rows if row["profile"] == profile
                             and row["V_over_2pi_MHz"] == interaction), key=lambda row: row["dt_s"])
            for ax, key in zip(axes, ("maximum_error_ratio", "final_fidelity_error")):
                ax.loglog([row["dt_s"] * 1e9 for row in subset], [row[key] for row in subset],
                          linestyle=style, marker=".", color=colors[index % len(colors)],
                          label=f"{profile}, {interaction:g} MHz")
    axes[0].axhline(1, color="0.4", linestyle=":", label="error threshold")
    axes[0].set_ylabel("Maximum error / tolerance")
    axes[1].set_ylabel("Final Bell fidelity error")
    axes[0].legend(fontsize=8, ncol=2)
    axes[-1].set_xlabel("Timestep (ns)")
    for ax in axes:
        ax.grid(which="both", alpha=0.15)
    fig.suptitle(case_label + "\nLeft-endpoint timestep convergence")
    save_figure(fig, figures / "sweep_convergence")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--interaction-mhz", type=float, nargs="+", required=True,
                        help="documented numerical samples of V/(2pi), in MHz")
    parser.add_argument("--milestone", action="store_true",
                        help="full 100 MHz benchmark plus cheaper family stability checks at a common timestep")
    parser.add_argument("--phase", type=float, required=True, help="constant laser phase in radians")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    parser.add_argument("--case-label", default="User-specified parameters",
                        help="parameter provenance label included in figures and report")
    args = parser.parse_args()
    if not np.all(np.isfinite(args.interaction_mhz)) or not np.isfinite(args.phase):
        parser.error("interaction and phase parameters must be finite")
    if len(set(args.interaction_mhz)) != len(args.interaction_mhz):
        parser.error("interaction values must be distinct")
    if args.milestone and (args.interaction_mhz != [0, 0.5, 1, 5, 100] or args.phase != 0):
        parser.error("--milestone requires --interaction-mhz 0 0.5 1 5 100 --phase 0")
    if args.milestone and args.case_label == "User-specified parameters":
        args.case_label = "Milestone: representative sampling of 0–100 MHz"
    figures, data = args.output_dir / "figures", args.output_dir / "data"
    figures.mkdir(parents=True, exist_ok=True)
    data.mkdir(parents=True, exist_ok=True)
    cases, all_rows, summaries, arrays = {}, [], [], {}
    production_steps = None
    benchmark_rows, benchmark_summaries = [], []
    if args.milestone:
        for profile in PROFILES:
            rows, summary, saved = run_case(profile, 100, args.phase, STEP_COUNTS)
            if summary["largest_tested_passing_dt_s"] is None:
                raise RuntimeError("Benchmark failed; refine STEP_COUNTS")
            benchmark_rows.extend(rows)
            benchmark_summaries.append(summary)
        production_steps = max(s["selected_steps"] for s in benchmark_summaries)
        # A common timestep is accepted only after every displayed case passes.
        for attempt in range(8):
            cases, all_rows, summaries = {}, [], []
            for profile in PROFILES:
                for interaction in args.interaction_mhz:
                    rows, summary, saved = run_case(profile, interaction, args.phase,
                                                   (2 * production_steps, production_steps),
                                                   family_check=interaction != 100)
                    cases[(profile, interaction)] = saved
                    all_rows.extend(rows)
                    summaries.append(summary)
            if all(s["selected_steps"] == production_steps and
                   s["largest_tested_passing_dt_s"] is not None for s in summaries):
                break
            production_steps *= 2
        else:
            raise RuntimeError("Family verification failed after eight refinements; figures withheld")
    for profile in PROFILES:
        for index, interaction in enumerate(args.interaction_mhz):
            if args.milestone:
                saved = cases[(profile, interaction)]
            else:
                rows, summary, saved = run_case(profile, interaction, args.phase, STEP_COUNTS)
                cases[(profile, interaction)] = saved
                all_rows.extend(rows)
                summaries.append(summary)
            arrays.update({f"{profile}_{index}_{key}": value for key, value in saved.items()})
    report = {
        "case_label": args.case_label,
        "source_interpretation": "Source specifies 0–100 MHz range, not individual red-curve values; milestone samples 0, 0.5, 1, 5, 100 MHz. Constant phase zero is a convention for phase-invariant probabilities from |00>.",
        "production_steps": production_steps,
        "production_dt_s": TMAX / production_steps if production_steps else None,
        "production_refinement_attempts": attempt + 1 if args.milestone else None,
        "family_criterion": "For each profile and V: P00 and P_Psi+ errors <=5e-5 against refined DOP853 at production and half timestep, and difference between those grids <=5e-5. The 100 MHz case additionally retains full scalar/trajectory criterion.",
        "benchmark_summaries": benchmark_summaries,
        "benchmark_convergence": benchmark_rows,
        "parameters": {"omega_rad_s": OMEGA, "delta_max_rad_s": DELTA_MAX, "tmax_s": TMAX,
                       "phase_rad": args.phase, "V_over_2pi_MHz": args.interaction_mhz,
                       "initial_state": "00", "basis_order": ["00", "0r", "r0", "rr"]},
        "method": "left-endpoint full matrix exponential, unchanged evolve_trotter",
        "criterion": "peak-scale half-unit fourth-significant-place trajectory errors; scalar final/peak four-digit rounding also identical",
        "reference": "DOP853 rtol=1e-13 atol=1e-15 max_step=T/4000; verified against rtol=1e-11 atol=1e-13 max_step=T/2000",
        "interpolation": "DOP853 dense output evaluated at each candidate's exact grid times; no interpolation of candidate probabilities",
        "peak_time": "candidate sampled maximum with dt resolution; reference dense maximum locally optimised after 40001-point scan; not part of acceptance criterion",
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "step_counts": STEP_COUNTS, "summaries": summaries, "convergence": all_rows,
    }
    (data / "sweep_validation.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    columns = ["profile", "V_over_2pi_MHz", "steps", "dt_s", "final_bell_fidelity",
               "final_fidelity_error", "maximum_P_00_error", "maximum_P_Psi_plus_error",
               "maximum_bell_fidelity_difference", "peak_fidelity_error", "maximum_error_ratio",
               "scalar_rounding_matches", "passes", "validation_scope",
               "halved_dt_maximum_plot_difference"]
    with (data / "sweep_convergence.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(benchmark_rows + [r for r in all_rows if r["validation_scope"] == "family plot stability"] if args.milestone else all_rows)
    np.savez_compressed(data / "figure_5_2.npz", **arrays,
                        interaction_mhz=args.interaction_mhz,
                        observable_labels=["P_00", "P_0r", "P_r0", "P_rr", "P_Psi_plus"])
    if any(summary["largest_tested_passing_dt_s"] is None for summary in summaries):
        raise RuntimeError("No passing timestep for at least one case; refine STEP_COUNTS. Data saved, figures withheld.")
    with plt.rc_context({"font.size": 10, "pdf.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False}):
        make_figures(cases, summaries, benchmark_rows + [r for r in all_rows if r["validation_scope"] == "family plot stability"] if args.milestone else all_rows, figures, args.case_label)
    print(json.dumps(summaries, indent=2), flush=True)


if __name__ == "__main__":
    main()

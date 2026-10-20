"""Plot recorded CSV data only. This script cannot run a benchmark."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def read_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def load_recorded_data(benchmark_dir, data_dir):
    """Reject mixed sessions or stale summaries without collecting any timings."""
    metadata = json.loads((benchmark_dir / "symmetry_benchmark_metadata.json").read_text())
    if not metadata.get("complete"):
        raise ValueError("Refusing to plot incomplete benchmark session")
    validation = json.loads((data_dir / "symmetry_validation.json").read_text())
    if (validation["benchmark_session_id"] != metadata["session_id"]
            or validation["exact_parameters"] != metadata["parameters"]):
        raise ValueError("Validation and benchmark data must belong to the same session")
    raw = read_rows(benchmark_dir / "symmetry_benchmark_raw.csv")
    rows = read_rows(benchmark_dir / "symmetry_benchmark_summary.csv")
    expected = {(n, "full") for n in metadata["full_sizes"]}
    expected.update((n, "reduced") for n in metadata["reduced_sizes"])
    groups = {key: [] for key in expected}
    for index, row in enumerate(raw):
        key = (int(row["N"]), row["method"])
        if (row["session_id"] != metadata["session_id"]
                or row["session_kind"] != metadata["session_kind"]
                or int(row["order"]) != index or key not in groups):
            raise ValueError("Raw timings contain inconsistent session or ordering metadata")
        times = [float(row[f"{phase}_s"]) for phase in ("construction", "evolution", "total")]
        if (not np.all(np.isfinite(times)) or min(times) < 0
                or not np.isclose(times[0] + times[1], times[2], rtol=1e-12, atol=1e-15)):
            raise ValueError("Raw timings contain invalid durations")
        groups[key].append(row)
    if len(rows) != len(expected) or {(int(r["N"]), r["method"]) for r in rows} != expected:
        raise ValueError("Summary sizes do not match the recorded session")
    for row in rows:
        n = int(row["N"])
        group = groups[(n, row["method"])]
        if (len(group) != metadata["repeats"] or int(row["repeats"]) != len(group)
                or {int(r["repeat"]) for r in group} != set(range(metadata["repeats"]))):
            raise ValueError("Raw timing repeats are missing or duplicated")
        for record in [row, *group]:
            if int(record["full_dimension"]) != 2**n or int(record["reduced_dimension"]) != n + 1:
                raise ValueError("Recorded dimensions do not match N")
        for phase in ("construction", "evolution", "total"):
            q25, median, q75 = np.percentile([float(r[f"{phase}_s"]) for r in group], [25, 50, 75])
            for suffix, value in (("median", median), ("q25", q25), ("q75", q75), ("iqr", q75-q25)):
                if not np.isclose(float(row[f"{phase}_{suffix}_s"]), value, rtol=1e-12, atol=1e-15):
                    raise ValueError("Summary statistics do not match raw timings")
        for field in ("hamiltonian_bytes", "trajectory_bytes", "estimated_working_bytes"):
            if int(row[field]) != max(int(r[field]) for r in group):
                raise ValueError("Summary memory does not match raw measurements")
    agreement = read_rows(data_dir / "symmetry_agreement.csv")
    if {int(r["N"]) for r in agreement} != set(metadata["full_sizes"]):
        raise ValueError("Agreement sizes do not match benchmark sizes")
    # The validation JSON duplicates these small metrics and ties them to a session.
    if len(agreement) != len(validation["exact_agreement"]):
        raise ValueError("Agreement CSV does not match validation metadata")
    for row, recorded in zip(agreement, validation["exact_agreement"]):
        if any(float(row[key]) != value for key, value in recorded.items()):
            raise ValueError("Agreement CSV does not match validation metadata")
    return metadata, rows, agreement, read_rows(data_dir / "symmetry_broken.csv")


def save(fig, directory, name):
    for suffix in ("png", "pdf"):
        fig.savefig(directory / f"{name}.{suffix}", dpi=180, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-dir", type=Path, default=ROOT / "results/data/symmetry_provisional")
    parser.add_argument("--data-dir", type=Path, default=ROOT / "results/data")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/figures")
    args = parser.parse_args()
    metadata, rows, agreement, broken = load_recorded_data(args.benchmark_dir, args.data_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    stamp = (f"PROVISIONAL · {metadata['machine_label']}"
             if metadata["session_kind"] == "provisional" else metadata["machine_label"])
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(7, 4), layout="constrained")
    n = np.arange(1, max(metadata["reduced_sizes"])+1)
    ax.semilogy(n, 2.**n, label=r"Full: $2^N$")
    ax.semilogy(n, n+1, label="Dicke: N+1")
    ax.set(xlabel="Number of atoms N", ylabel="Hilbert-space dimension", title="S1 · Hilbert-space reduction")
    ax.grid(alpha=.25); ax.legend()
    save(fig, args.output_dir, "S1_dimensions")

    for phases, name in [(["total"], "S2_runtime"), (["construction", "evolution"], "S2_runtime_components")]:
        fig, axes = plt.subplots(1, len(phases), figsize=(7*len(phases), 4), layout="constrained", squeeze=False)
        for ax, phase in zip(axes.flat, phases):
            for method, label in (("full", "Full computational basis"), ("reduced", "Dicke basis")):
                group = [r for r in rows if r["method"] == method]
                ns = [int(r["N"]) for r in group]
                med = np.array([float(r[f"{phase}_median_s"]) for r in group])
                low = np.array([float(r[f"{phase}_q25_s"]) for r in group])
                high = np.array([float(r[f"{phase}_q75_s"]) for r in group])
                ax.errorbar(ns, med, yerr=[med-low, high-med], fmt="o-", capsize=3, label=label)
            ax.set_yscale("log")
            ax.set_xscale("log", base=2)
            ax.set_xticks([2, 4, 8, 12, 16, 32, 64], ["2", "4", "8", "12", "16", "32", "64"])
            ax.set(xlabel="Number of atoms N", ylabel=f"{phase.capitalize()} time (s)")
            ax.grid(alpha=.25); ax.legend()
        fig.suptitle(f"S2 · Median and IQR, {metadata['repeats']} repeats · {stamp}")
        save(fig, args.output_dir, name)

    fig, ax = plt.subplots(figsize=(7, 4), layout="constrained")
    for phase in ("total", "construction", "evolution"):
        sizes, speed = [], []
        for n in metadata["full_sizes"]:
            lookup = {r["method"]: r for r in rows if int(r["N"]) == n}
            sizes.append(n)
            speed.append(float(lookup["full"][f"{phase}_median_s"]) / float(lookup["reduced"][f"{phase}_median_s"]))
        ax.semilogy(sizes, speed, "o-", label=phase.capitalize())
    ax.axhline(1, linestyle=":", color="0.5")
    ax.set(xlabel="Number of atoms N", ylabel="Median full time / median reduced time", title=f"S3 · Speed-up · {stamp}")
    ax.grid(alpha=.25); ax.legend()
    save(fig, args.output_dir, "S3_speedup")

    fig, ax = plt.subplots(figsize=(7, 4), layout="constrained")
    ns = [int(r["N"]) for r in agreement]
    for key, label in (("raw_infidelity", "1 − minimum fidelity (raw)"), ("max_dicke_population_error", "Maximum Dicke population error")):
        values = np.array([float(r[key]) for r in agreement])
        ax.plot(ns, values, "o-", label=label)
    ax.set_yscale("symlog", linthresh=1e-15)
    ax.axhline(0, color="0.5", linestyle=":", label="Zero; signed roundoff retained")
    ax.set(xlabel="Number of atoms N", ylabel="Numerical discrepancy", title="S4 · Exact symmetry: full versus Dicke")
    ax.grid(alpha=.25); ax.legend(fontsize=9)
    save(fig, args.output_dir, "S4_agreement")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    for n in sorted({int(r["N"]) for r in broken}):
        group = [r for r in broken if int(r["N"]) == n]
        x = [float(r["alpha"]) for r in group]
        for ax, key in zip(axes, ("infidelity", "max_symmetric_leakage")):
            ax.plot(x, [float(r[key]) for r in group], "o-", label=f"N={n}")
    for ax, label in zip(axes, ("1 − minimum fidelity", "Maximum leakage from Dicke space")):
        ax.set(xlabel="Nonuniformity interpolation α (1 = physical chain)", ylabel=label)
        ax.grid(alpha=.25)
    axes[0].legend()
    fig.suptitle("S5 · Breakdown of the uniform-interaction approximation")
    save(fig, args.output_dir, "S5_broken_symmetry")

    fig, ax = plt.subplots(figsize=(7, 4), layout="constrained")
    for method in ("full", "reduced"):
        group = [r for r in rows if r["method"] == method]
        ax.semilogy([int(r["N"]) for r in group],
                    [(float(r["hamiltonian_bytes"])+float(r["trajectory_bytes"]))/2**20 for r in group],
                    "o-", label=f"{method}: CSR Hamiltonian + trajectory")
    ax.set(xlabel="Number of atoms N", ylabel="Array payload (MiB; excludes workspace)", title="S6 · Stored memory for the benchmark grid")
    ax.grid(alpha=.25); ax.legend(fontsize=9)
    save(fig, args.output_dir, "S6_memory")
    print(f"Read recorded data from {args.benchmark_dir}; wrote S1–S6 to {args.output_dir}")


if __name__ == "__main__":
    main()

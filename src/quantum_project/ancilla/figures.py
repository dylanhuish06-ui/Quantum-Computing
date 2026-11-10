"""Publication figures and explicit cross-model error taxonomy, from raw data."""
import csv
import json
from pathlib import Path
import uuid

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .run import digest, environment, sha256

COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]
STYLES = ["-", "--", "-.", ":"]


def valid_products(directory, spec):
    try:
        manifest = json.loads((directory/"manifest.json").read_text())
        if manifest["specification"] != spec or manifest["input_hash"] != digest(spec):
            return False
        hashes = manifest["artifact_sha256"]
        required = {"report.md", "errors.csv", "figure_data.npz", "comparison.json"}
        required |= {name+ext for name in ("accuracy", "protection", "large_n", "adaptive") for ext in (".pdf", ".png")}
        if not required.issubset(hashes):
            return False
        if {p.name for p in directory.iterdir() if p.is_file()} != set(hashes)|{"manifest.json"}:
            return False
        return all(Path(name).name == name and sha256(directory/name) == checksum
                   for name, checksum in hashes.items())
    except (OSError, ValueError, KeyError, TypeError):
        return False


def error_table(reports):
    rows = []
    def add(category, N, resource, metric, value, meaning):
        rows.append(dict(category=category, N=N, resource=resource, metric=metric,
                         value=value, interpretation=meaning))
    for r in reports[4]["rows"]:
        category = "exact symmetry" if r["m"] == 0 else "exact hybrid structured asymmetry"
        add(category, r["N"], f"m={r['m']}", "complex_trajectory_l2", r["complex_trajectory_l2"],
            "numerical difference from independent full-space evolution; no phase alignment")
    for r in reports[2]["rows"]:
        for key in ("conditional_fidelity", "success_probability", "successful_target_probability"):
            add("postselected protected dynamics", r["N"], f"checks={r['checks']}", key, r[key],
                "physical protocol outcome relative to projected target; not numerical error")
        add("postselected protected dynamics", r["N"], f"checks={r['checks']}",
            "ode_refinement_l2", r["ode_refinement_l2"], "numerical integration refinement at fixed checks")
    for r in reports[3]["rows"]:
        add("coherent gap protection", r["N"], f"G/Omega_ref={r['G_over_Omega_ref']}",
            "final_leakage", r["final_leakage"], "physical departure from symmetry; not numerical error")
    for r in reports[5]["rows"]:
        for key in ("estimated_fine_time_error_l2", "action_refinement_l2", "output_grid_l2"):
            add("exact hybrid structured asymmetry", r["N"], "m=1; large N", key, r[key],
                "deterministic numerical diagnostic; time error is asymptotic Richardson estimate")
    for r in reports[5]["adaptive"]:
        if rows and any(x["category"] == "approximate adaptive dynamics" and x["resource"] == f"m={r['m']}" for x in rows):
            continue
        for key in ("state_error_bound", "actual_state_error_l2"):
            add("approximate adaptive dynamics", 5, f"m={r['m']}", key, r[key],
                "static projected-model error; Duhamel bound excludes roundoff; full-space diagnostic")
    return rows


def render(reports, directories, output):
    plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11,
                         "legend.fontsize": 9, "pdf.fonttype": 42, "ps.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.facecolor": "white"})
    def save(fig, name):
        for ext in ("pdf", "png"):
            fig.savefig(output/f"{name}.{ext}", dpi=220, bbox_inches="tight")
        plt.close(fig)
    data = {}
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), layout="constrained")
    for m in range(4):
        rows = [r for r in reports[4]["rows"] if r["m"] == m]
        label = "Symmetric (m=0)" if m == 0 else f"Hybrid (m={m})"
        axes[0].semilogy([r["N"] for r in rows], np.maximum([r["complex_trajectory_l2"] for r in rows], 1e-16),
                         STYLES[m], color=COLORS[m], marker="o", markersize=4, label=label)
    axes[0].set(xlabel="Total atom count N", ylabel="Maximum complex-state L2 difference",
                title="a  Exact reductions vs independent full space")
    axes[0].legend()
    axes[0].text(.03, .03, r"Zero residuals displayed at $10^{-16}$", transform=axes[0].transAxes, fontsize=8)
    rows = reports[5]["rows"]
    for i, r in enumerate(rows):
        axes[1].loglog([.01, .005], [r["coarse_medium_l2"], r["medium_fine_l2"]],
                       STYLES[i], color=COLORS[i], marker="o", label=f"N={r['N']}; order {r['midpoint_order']:.3f}")
    axes[1].set(xlabel=r"Coarser midpoint step $\Delta\tau$", ylabel="Difference between successive refinements",
                title="b  Large-N time-step convergence")
    axes[1].legend()
    for ax in axes:
        ax.grid(alpha=.2)
    save(fig, "accuracy")

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), layout="constrained")
    zeno = [r for r in reports[2]["rows"] if r["N"] == 6]
    for i, (key, label) in enumerate((("conditional_fidelity", "Conditional target fidelity"),
                                    ("success_probability", "All-check success probability"),
                                    ("successful_target_probability", "Successful target probability"))):
        values = [r[key] for r in zeno]
        axes[0].plot([r["checks"] for r in zeno], values, STYLES[i], color=COLORS[i],
                     marker=("o", "s", "^")[i], markerfacecolor="none", label=label)
        data["zeno_"+key] = values
    axes[0].set(xscale="log", xlabel="Number of ideal syndrome checks", ylabel="Probability",
                ylim=(min(r["successful_target_probability"] for r in zeno)-.015, 1.003),
                title=r"a  Selective protection: N=6, $\tau=1.6$")
    axes[0].legend(loc="lower right")
    gap = [r for r in reports[3]["rows"] if r["N"] == 6]
    axes[1].semilogy([r["G_over_Omega_ref"] for r in gap], [r["final_leakage"] for r in gap],
                     "o-", color=COLORS[0], label="Coherent spin-gap protection")
    axes[1].axhline(zeno[0]["bare_leakage"], color=COLORS[1], ls="--", label="Unprotected baseline")
    axes[1].set(xlabel=r"Penalty strength $G/\Omega_{\rm ref}$", ylabel="Final symmetry leakage",
                title="b  Same disorder, duration and initial state")
    axes[1].legend()
    for ax in axes:
        ax.grid(alpha=.2)
    save(fig, "protection")

    with (directories[5]/"observables.csv").open() as stream:
        obs = list(csv.DictReader(stream))
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), layout="constrained")
    for i, N in enumerate((64, 256, 1024)):
        rows = [r for r in obs if int(r["N"]) == N]
        t = np.array([float(r["time_seconds"])*1e9 for r in rows])
        for ax, key in zip(axes, ("bulk_excitation_fraction", "defect_population")):
            values = np.array([float(r[key]) for r in rows])
            ax.plot(t, values, STYLES[i], color=COLORS[i], linewidth=2, label=f"N={N}")
            data[f"N{N}_"+key] = values
        data[f"N{N}_time_ns"] = t
    for ax, title, ylabel in zip(axes, ("a  Symmetric bulk", "b  Distinguished physical atom"),
                                ("Excitation fraction per bulk atom", "Excitation probability")):
        ax.set(xlabel="Time (ns)", ylabel=ylabel, title=title)
        ax.grid(alpha=.2)
        ax.legend()
    save(fig, "large_n")

    adaptive = {r["m"]: r for r in reports[5]["adaptive"]}
    rows = list(adaptive.values())
    fig, ax = plt.subplots(figsize=(7.2, 4.2), layout="constrained")
    for i, (key, label) in enumerate((("state_error_bound", "Duhamel upper bound"),
                                    ("actual_state_error_l2", "Measured model error"))):
        ax.plot([r["dimension"] for r in rows], [r[key] for r in rows], STYLES[i],
                color=COLORS[i], marker="o", label=label)
        data["adaptive_"+key] = [r[key] for r in rows]
    ax.axhline(1e-4, ls=":", color=COLORS[2], label="Required error budget: 0.0001")
    ax.set(xlabel="Retained state dimension", ylabel="Complex-state L2 error / bound",
           title="Generic disorder: tight budget requires full space (N=5)")
    ax.set_xticks([r["dimension"] for r in rows])
    ax.legend()
    ax.grid(alpha=.2)
    save(fig, "adaptive")
    np.savez_compressed(output/"figure_data.npz", **data)


def build_figures(directories, root, *, resume):
    reports = {k: json.loads((v/"validation.json").read_text()) for k, v in directories.items()}
    spec = {"schema": 1, "source_sha256": sha256(Path(__file__)), "environment": environment(),
            "experiment_manifests": {str(k): sha256(v/"manifest.json") for k, v in directories.items()}}
    directory = root/("figures-"+digest(spec)[:16])
    if directory.exists():
        if resume and valid_products(directory, spec):
            print("Figures: verified cache, skipped", flush=True)
            return directory
        if not resume:
            raise ValueError("existing figure outputs preserved; use --resume")
        directory.rename(root/(directory.name+".invalid-"+uuid.uuid4().hex[:8]))
    pending = root/(directory.name+".pending-"+uuid.uuid4().hex[:8])
    pending.mkdir()
    render(reports, directories, pending)
    errors = error_table(reports)
    with (pending/"errors.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(errors[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(errors)
    comparison = {"schema": 1, "input_specification": spec, "errors": errors,
                  "experiment_parameters": {str(k): r["parameters"] for k, r in reports.items()},
                  "experiment_environment": {str(k): r["environment"] for k, r in reports.items()},
                  "experiment_provenance": {str(k): json.loads((directories[k]/"manifest.json").read_text()) for k in reports},
                  "uncertainty": "deterministic solver/model diagnostics; no sampling error bars",
                  "benchmark": "same-session environment; no runtime speedup claims"}
    (pending/"comparison.json").write_text(json.dumps(comparison, indent=2, allow_nan=False)+"\n")
    symmetric = max(r["complex_trajectory_l2"] for r in reports[4]["rows"] if r["m"] == 0)
    hybrid = max(r["complex_trajectory_l2"] for r in reports[4]["rows"] if r["m"] > 0)
    numerical_zeno = max(r["ode_refinement_l2"] for r in reports[2]["rows"])
    largest = reports[5]["rows"][-1]
    text = f"""# Phase 4 comparison and reproducibility

All experiments were recomputed in this session environment. Raw CSV/NPZ data,
parameter recipes, software versions, local Git provenance and content hashes
are retained in the experiment directories linked by `latest.json`.

| Class | Validation / result | Interpretation |
| --- | --- | --- |
| Exact symmetry | Maximum full-reference L2 difference {symmetric:.3e} | Bulk symmetry and symmetric initial bulk required |
| Exact hybrid structured asymmetry | Maximum small-N reference difference {hybrid:.3e} | Distinguished physical atoms; invariant bulk |
| Postselected protected dynamics | Maximum ODE refinement difference {numerical_zeno:.3e} | Physical target is projected dynamics; success probability is reported separately |
| Approximate adaptive dynamics | N=5 tight budget reaches dimension 32 | Static full-space diagnostic; generic compression can fail |

The large-N model uses N=64,256,1024, one distinguished physical atom and
tau=Omega_ref*t with Omega_ref=2*pi*1 MHz. At N=1024 the finest time-step error
estimate is {largest['estimated_fine_time_error_l2']:.3e}; norm drift is
{largest['norm_drift']:.3e}. The time error is an asymptotic Richardson
estimate, not a rigorous error bar. Exponential action, time step and output
sampling are refined independently. Initial states and complete coefficient
profiles are specified by the content-hashed validation recipes.

`errors.csv` separates numerical residuals, physical leakage, postselection
probabilities and approximation error. These quantities must not be ranked as
if they measured the same task. They are deterministic; no statistical error
bars or experimental noise estimates are inferred. `comparison.json` retains
the full machine-readable table and provenance.

`accuracy.pdf/png` compares exact reductions with independent full-space
evolution and shows large-N second-order convergence. Exactly zero small-N
residuals are displayed at a labelled 1e-16 plotting floor; raw values are
retained unchanged in the error table. `protection.pdf/png`
compares ideal selective checks and coherent spin-gap protection at the same
N=6 disorder, duration and initial state. Conditional fidelity, all-check
success and successful-target probability are different observables.
`large_n.pdf/png` shows bulk and distinguished-site populations; almost
overlapping curves are a measured feature of this short sweep.
`adaptive.pdf/png` shows the conservative static Duhamel residual bound and
actual model error; neither bound nor error must decrease monotonically with
each promoted atom. Reaching dimension 32 eliminates compression at N=5.

The ideal syndrome oracle is nonlocal; its gate compilation depth and hardware
cost remain unresolved. The coherent total-spin penalty is not hardware-free.
No local hardware implementation, generic exact polynomial simulator or
speedup benchmark is claimed. CPU timings are provisional runtime records.

Repeat `PYTHONPATH=src python -m quantum_project.ancilla.run --resume` from the
repository. Only matching parameter/source/environment hashes and intact
artifact checksums are skipped. Invalid outputs are preserved with an invalid
suffix and rebuilt; interrupted pending folders remain available for inspection.
This is an experiment runner, not an autonomous code-generation command.
"""
    (pending/"report.md").write_text(text)
    manifest = {"input_hash": digest(spec), "specification": spec,
                "artifact_sha256": {p.name: sha256(p) for p in pending.iterdir() if p.is_file()}}
    (pending/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    if not valid_products(pending, spec):
        raise ValueError("figure manifest validation failed")
    pending.rename(directory)
    print("Figures: generated and verified", flush=True)
    return directory

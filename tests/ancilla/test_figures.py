import json

from quantum_project.ancilla.figures import valid_products, error_table
from quantum_project.ancilla.run import digest, sha256


def test_figures_require_every_expected_artifact_and_checksum(tmp_path):
    spec = {"parameters": "fixed"}
    names = ["report.md", "errors.csv", "figure_data.npz", "comparison.json"]
    names += [name+ext for name in ("accuracy", "protection", "large_n", "adaptive") for ext in (".pdf", ".png")]
    for name in names:
        (tmp_path/name).write_bytes(b"test content")
    manifest = {"specification": spec, "input_hash": digest(spec),
                "artifact_sha256": {name: sha256(tmp_path/name) for name in names}}
    (tmp_path/"manifest.json").write_text(json.dumps(manifest))
    assert valid_products(tmp_path, spec)
    assert not valid_products(tmp_path, {"parameters": "changed"})
    (tmp_path/"protection.pdf").write_bytes(b"altered")
    assert not valid_products(tmp_path, spec)
    (tmp_path/"protection.pdf").unlink()
    assert not valid_products(tmp_path, spec)


def test_error_table_distinguishes_physical_and_numerical_quantities():
    reports = {
        2: {"rows": [{"N": 2, "checks": 8, "conditional_fidelity": .99,
                       "success_probability": .8, "successful_target_probability": .792,
                       "ode_refinement_l2": 1e-13}]},
        3: {"rows": [{"N": 2, "G_over_Omega_ref": 1, "final_leakage": .01}]},
        4: {"rows": [{"N": 2, "m": m, "complex_trajectory_l2": 1e-12} for m in (0, 1)]},
        5: {"rows": [{"N": 64, "estimated_fine_time_error_l2": 1e-5,
                       "action_refinement_l2": 1e-14, "output_grid_l2": 1e-15}],
            "adaptive": [{"m": 0, "state_error_bound": .4, "actual_state_error_l2": .1}]}}
    rows = error_table(reports)
    assert len({r["category"] for r in rows}) == 5
    physical = [r for r in rows if r["metric"] == "conditional_fidelity"][0]
    assert "not numerical error" in physical["interpretation"]
    assert all("numerical" in r["interpretation"] for r in rows if r["metric"] == "complex_trajectory_l2")

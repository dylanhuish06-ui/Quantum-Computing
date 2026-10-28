"""Read-only audit of the selected Phase 4 checkpoint; run with PYTHONPATH=src."""
import json
from pathlib import Path
import re

import numpy as np
from quantum_project.ancilla.run import ROOT, specification, valid_cache
from quantum_project.ancilla.figures import valid_products


def audit():
    root = ROOT/"results/ancilla/runs"
    latest = json.loads((root/"latest.json").read_text())
    checked, arrays = [], 0
    for gate, name in latest["gates"].items():
        directory = root/name
        assert valid_cache(directory, specification(int(gate))), name
        for path in directory.glob("*.npz"):
            with np.load(path, allow_pickle=False) as data:
                for key in data.files:
                    assert np.all(np.isfinite(data[key])), (path, key)
                    arrays += 1
        checked.append(name)
    directory = root/latest["figures"]
    manifest = json.loads((directory/"manifest.json").read_text())
    assert valid_products(directory, manifest["specification"])
    for gate, name in latest["gates"].items():
        from quantum_project.ancilla.run import sha256
        assert manifest["specification"]["experiment_manifests"][gate] == sha256(root/name/"manifest.json")
    for name in ("accuracy", "adaptive", "large_n", "protection"):
        assert (directory/(name+".pdf")).read_bytes().startswith(b"%PDF-")
        assert (directory/(name+".png")).read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    links = []
    for path in (ROOT/"docs/ancilla_phase4").glob("*.md"):
        for target in re.findall(r"\]\(([^)]+)\)", path.read_text()):
            if "://" in target or target.startswith("#"):
                continue
            resolved = path.parent/target.split("#")[0]
            assert resolved.exists(), (path, target)
            links.append(str(resolved.relative_to(ROOT)))
    return {"passed": True, "verified_gates": checked,
            "finite_trajectory_arrays": arrays, "figure_pairs": 4,
            "documentation_links": links,
            "checks": ["scientific source/environment signatures", "artifact inventories and SHA256",
                       "cross-manifest figure provenance", "finite raw trajectories", "PDF/PNG signatures", "relative documentation links"]}


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))

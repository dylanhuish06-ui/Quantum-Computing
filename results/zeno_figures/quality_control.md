# Final quality-control record

Nine numbered scientific figures, ten export sets, and PDF/PNG contact sheets are complete. No poster or presentation was created.

## Tests and independent reproduction

- Original working-tree tests before additions: **364 passed**.
- Full working tree after additions: **381 passed**, no failures, errors or skips (15.45 s).
- Clean copy of tracked HEAD plus only this figure-suite addition: **305 passed**, no failures, errors or skips (11.94 s). The additional 76 working-tree tests belong to the pre-existing untracked Zeno extension.
- Regeneration in that clean copy reproduced **47 files byte for byte**, including the numerical archive, newly computed dataset, every export, manifest, guide, requirements and reproducibility report. The clean copy used the same scientific virtual environment but contained no untracked original Zeno modules, tests or data. Manual QA files were written after this comparison and are not auto-generated artifacts.
- All 90 original files fingerprinted before work remain byte-identical. The original index was empty. No frozen milestone or original research dataset was changed.
- The manual GitHub Actions workflow has not been run on GitHub. Its generation and test commands were exercised locally; no push was performed.

## Numerical checks

The complete saved archive was independently checked, including 24,192 complex density matrices. Maximum saved trace error: 8.340e-13; minimum density eigenvalue: -2.442e-15 (roundoff). The largest reference discrepancy was 1.431e-11, below the explicit 3e-11 tolerance for accumulated high-count propagation. Hamiltonian Hermiticity, exact projector definitions, complete non-selective channel, M=0/M=1, outcome-order invariance, short-time quadratic coefficient, selective/non-selective ordering and large-M leading loss all passed.

The Hamiltonian gives tau_Z=225.079079039 ns and a=(T/tau_Z)^2=4.934802200545. Independent tests additionally enumerate all outcome histories, recover the V=0 analytical all-success probability, check fourth-order short-time residuals, and verify dynamics inside the rank-three protected sector. RK4 convergence approaches fourth order; the maximum exact-exponential subdivision discrepancy is 4.778e-14. No probability clipping, artificial density repair, error bars or experimental claims were introduced.

## Visual inspection actually performed

All ten PDFs were rendered with **PyMuPDF 1.28.2 at 130 DPI**, and the resulting images were opened and inspected. The final contact-sheet PDF was rendered and inspected at 100 DPI. Figures 01, 05, 07 and 08 were inspected again following layout/rendering corrections; the single-column Figure 03 was checked separately. Review covered legibility, legends, mathematical symbols, scales, axis units, parameter annotations, label bounds and complete contact-sheet coverage.

Each numbered PDF has one page, embedded fonts and **zero embedded raster images**. Each SVG parses successfully, retains text and contains no image elements. Numbered PNGs have 400-DPI metadata: 3360 pixels across for wide figures and 1680 for the single-column variant. The contact PNG is 3600 by 4800 pixels at 300 DPI. Every export hash and individual raster dimension appears in the manifest; PDF structure checks and test summaries appear in `quality_control.json`.

Initial export attempts exposed a math-text escape error, a missing non-math glyph, rasterised colour bars, a clipped heatmap label and cramped narrow-figure numbering. These were corrected. PDF inspection then exposed box-edge clipping, colour-mesh seams and an uneven convergence-panel heading; these were also corrected. The final exports passed the label-bound check and visual review. Two background verification attempts were interrupted with termination signals, without completed test results; sequential reruns completed successfully. No failed test is outstanding. Git whitespace checking initially flagged Matplotlib’s canonical SVG path-line spaces. A figure-directory `.gitattributes` now exempts only generated SVG whitespace and marks those files as generated; their deterministic bytes are preserved. The resulting staged diff check passes.

## Scope limits

Optional Figure 10 was omitted because the inspected symmetry worktree contains validated reduction/leakage data but no validated projection-protected large-N frequency study. All graphics describe ideal instantaneous projections in a static two-atom model. The contact sheet is only a selection aid. PDF is the preferred portable vector export; SVG live text depends on DejaVu/STIX font availability in the target editor.

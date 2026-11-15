# Missed-candidate audit

Reanalyse saved 94-candidate rankings against pinned source MLC readings and
raw annotation pointers. No training or model inference.

[Findings, scope and reproduction](../../docs/missed-candidate-audit.md).
[Versioned results](../../data/studies/results/missed_candidate_audit.v1/).

`build_audit.py` uses the Python 3.11+ standard library and verifies every input
against `input_manifest.json`. `verify_audit.py` independently checks source cells,
top-two decisions and rank intervals without importing the builder.

Downloads and locally generated raw annotation extracts belong under
`data/raw/` and `artifacts/`, respectively. Use the reproduction commands in the
report; do not commit those source files or local workbook intermediates.

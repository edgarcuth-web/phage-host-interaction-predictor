"""Compare source-defined discovery routes on a locked, observed phage panel.

No GPU, external model, provider API, sequencing service or wet-lab execution.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from phage_cocktail.retrieval_rescue import run_discovery_audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True,
                        help="Local v1 JSON manifest with source-hashed CSV inputs")
    parser.add_argument("--outdir", type=Path, required=True,
                        help="New, non-existing directory for immutable audit results")
    args = parser.parse_args()
    report = run_discovery_audit(args.manifest, args.outdir)
    print(f"Study: {report['study_id']}")
    print(f"Endpoint: {report['endpoint']}")
    print(f"Complete hosts: {report['method_summaries']['evidence']['complete_panel_hosts']}")
    print(f"Learned-only observed positives at K={report['primary_k']}: "
          f"{report['primary_k_rescue']['learned_only_positive_instances']}")
    print(f"Union-only observed positives at K={report['primary_k']}: "
          f"{report['primary_k_rescue']['union_only_positive_instances']}")
    print(f"Report: {args.outdir / 'report.json'}")


if __name__ == "__main__":
    main()

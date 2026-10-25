"""Reproduce the published retrieval metrics from locked precomputed matrices."""

import argparse
from pathlib import Path

from phage_cocktail.evo2_reference import run_reference


def main():
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true", help="Download missing locked upstream files")
    parser.add_argument("--source-root", type=Path, default=repo / "external/Evo2-Host-Prediction")
    parser.add_argument("--outdir", type=Path, default=repo / "artifacts/evo2_retrieval_reference_v1")
    args = parser.parse_args()
    report = run_reference(repo, args.source_root, args.outdir, fetch=args.fetch)
    print(f"Report: {args.outdir.resolve() / 'report.json'}")
    if not report["matches_reference"]:
        raise SystemExit("Reference metric mismatch; inspect report.json")
    print("All six methods match the rounded upstream reference metrics.")


if __name__ == "__main__":
    main()

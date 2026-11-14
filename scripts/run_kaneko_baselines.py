"""Fit or evaluate the frozen Kaneko known-library research baselines."""
import argparse
from pathlib import Path

from phage_cocktail.kaneko_baseline import evaluate_study, fit_study


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="stage", required=True)
    for name in ("fit", "evaluate"):
        stage = sub.add_parser(name)
        stage.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
        stage.add_argument("--grades", type=Path, required=True)
        stage.add_argument("--out-dir", type=Path, required=True)
    fit = sub.choices["fit"]
    fit.add_argument("--sequence-root", type=Path, required=True)
    fit.add_argument("--implementation-commit", required=True)
    evaluation = sub.choices["evaluate"]
    evaluation.add_argument("--fit-dir", type=Path, required=True)
    evaluation.add_argument("--frozen-fit-commit", required=True)
    args = parser.parse_args()
    if args.stage == "fit":
        fit_study(args.root, args.sequence_root, args.grades, args.out_dir, args.implementation_commit)
    else:
        evaluate_study(args.root, args.fit_dir, args.grades, args.out_dir, args.frozen_fit_commit)


if __name__ == "__main__":
    main()

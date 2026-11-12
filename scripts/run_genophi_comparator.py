"""Run or resume the prespecified GenoPHI comparator on the frozen host split."""
import argparse
from pathlib import Path

from phage_cocktail.genophi_comparator import run_comparator


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    parser.add_argument("--proteomes", type=Path, default=Path("data/derived/genophi_proteomes"))
    parser.add_argument("--proteome-receipt", type=Path, default=Path("data/receipts/genophi_proteomes.v1.json"))
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--stage", choices=["prepare", "cluster", "fit", "all"], default="all")
    parser.add_argument("--scenario", choices=["grouped_host", "phylogroup:B2"])
    parser.add_argument("--repeat", type=int, choices=[0, 1, 2])
    parser.add_argument("--protocol", type=Path, help="Frozen v1 or registered seed-43/44 comparator protocol")
    parser.add_argument("--feature-cache", type=Path, help="Hash-verified v1 feature-stage source for sensitivities")
    args = parser.parse_args()
    run_comparator(Path(__file__).resolve().parents[1], args.data_root, args.proteomes,
                   args.proteome_receipt, args.outdir, stage=args.stage, scenario=args.scenario, repeat=args.repeat,
                   protocol=args.protocol, feature_cache=args.feature_cache)


if __name__ == "__main__":
    main()

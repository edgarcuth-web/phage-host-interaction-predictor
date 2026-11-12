"""Compare host genome composition and annotation features on a fixed phage panel."""
import argparse
from pathlib import Path

from phage_cocktail.retrieval_study import run_retrieval_study


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--fetch", action="store_true")
    parser.add_argument("--protocol", type=Path, help="Registered frozen seed-sensitivity protocol")
    parser.add_argument("--scenario", choices=["grouped_host"], help="Execute only the grouped-host sensitivity")
    args = parser.parse_args()
    run_retrieval_study(Path(__file__).resolve().parents[1], args.data_root,
                        args.outdir, fetch=args.fetch, genome=True, protocol=args.protocol, scenario=args.scenario)


if __name__ == "__main__":
    main()

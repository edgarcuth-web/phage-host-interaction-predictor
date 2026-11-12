"""Freeze the registered seed-43/44 manifests and strictly inherited protocols; no fitting."""
import argparse
from pathlib import Path

from phage_cocktail.genophi_seed_protocols import freeze_seed_inputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    parser.add_argument("--feature-cache", type=Path, required=True)
    parser.add_argument("--request-audit", type=Path, required=True)
    args = parser.parse_args()
    result = freeze_seed_inputs(Path(__file__).resolve().parents[1], args.data_root,
                                args.feature_cache, args.request_audit)
    print(result["seed_splits"])


if __name__ == "__main__":
    main()

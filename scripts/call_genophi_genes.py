"""Call code-11 proteins from the verified frozen GenoPHI nucleotide cohort."""
import argparse
from pathlib import Path

from phage_cocktail.genophi_proteomes import prepare_proteomes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    parser.add_argument("--proteomes", type=Path, default=Path("data/derived/genophi_proteomes"))
    parser.add_argument("--out", type=Path, required=True, help="New immutable proteome receipt")
    parser.add_argument("--source-receipt", type=Path, help="Also record all verified source hashes")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    result = prepare_proteomes(
        Path(__file__).resolve().parents[1], args.data_root, args.proteomes, args.out,
        workers=args.workers, source_receipt=args.source_receipt,
    )
    print(result["counts"])
    print(f"Wall time: {result['wall_seconds']:.1f} s")


if __name__ == "__main__":
    main()

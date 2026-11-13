#!/usr/bin/env python3
"""Check the frozen retrieval study inputs; never fit a model or score a pair."""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from phage_cocktail.retrieval_protocol import check_protocol


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--annotation-data-root", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    try:
        if args.out is not None and args.out.exists():
            raise ValueError("output must be new")
        result = check_protocol(args.root, annotation_data_root=args.annotation_data_root)
        payload = json.dumps(result, indent=2) + "\n"
        if args.out is not None:
            args.out.parent.mkdir(parents=True, exist_ok=True)
            with args.out.open("x") as stream:
                stream.write(payload)
        print(payload, end="")
    except (ValueError, OSError, KeyError, TypeError) as error:
        print(f"preflight failed: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Audit potential external phage–host evaluation row ledgers; no data is fetched."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from phage_cocktail.external_overlap import REGIMES, audit_overlap, read_intake


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--development", type=Path, required=True, help="CSV development row ledger")
    parser.add_argument("--candidate", type=Path, required=True, help="CSV proposed external row ledger")
    parser.add_argument("--regime", choices=sorted(REGIMES), required=True)
    parser.add_argument("--out", type=Path, required=True, help="New JSON report path (will not overwrite)")
    args = parser.parse_args()
    report = audit_overlap(read_intake(args.development), read_intake(args.candidate), regime=args.regime)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x", encoding="utf-8") as file:
        json.dump(report, file, indent=2, sort_keys=True)
        file.write("\n")
    print(f"{report['status']}: {args.out}")
    return 2 if report["status"] == "blocked" else 0


if __name__ == "__main__":
    raise SystemExit(main())

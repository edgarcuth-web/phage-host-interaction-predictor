"""Run or resume the frozen three-host hosted-API engineering pilot."""
import argparse
import getpass
import os
from pathlib import Path

from phage_cocktail.nvidia_host_pilot import run_pilot


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    parser.add_argument("--feasibility-cache", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--verify-only", action="store_true", help="Validate a completed cache without a key or network")
    args = parser.parse_args()
    if args.verify_only and not (args.outdir / "report.json").exists():
        raise ValueError("verify-only requires a completed pilot report")
    key = None if args.verify_only else (os.environ.get("NVIDIA_API_KEY") or getpass.getpass("NVIDIA API key (hidden): "))
    try:
        report = run_pilot(Path(__file__).resolve().parents[1], args.data_root, args.outdir,
                           key, args.feasibility_cache, resume=args.resume or args.verify_only)
    finally:
        key = None
    print({k: report[k] for k in ("status", "hosts", "host_vector_shape", "new_post_attempts", "new_response_bytes")})


if __name__ == "__main__":
    main()

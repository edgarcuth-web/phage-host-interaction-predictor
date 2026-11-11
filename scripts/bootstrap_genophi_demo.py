"""Fetch only GenoPHI's pinned test-data directory using sparse checkout."""

from __future__ import annotations

from pathlib import Path
import shutil
import subprocess

from phage_cocktail.genophi_reference import GENOPHI_COMMIT, GENOPHI_REPOSITORY

TARGET = Path("external/GenoPHI")
REQUIRED_RELATIVE = (
    Path("data/test_data/strain_AAs"),
    Path("data/test_data/phage_AAs"),
    Path("data/test_data/ecoli_test_interaction_matrix.csv"),
)


def worktree_is_populated() -> bool:
    """True only when every required demo path is actually checked out."""
    return all((TARGET / relative).exists() for relative in REQUIRED_RELATIVE)


def run(*args: str) -> None:
    subprocess.run(args, check=True)


def current_commit() -> str | None:
    if not (TARGET / ".git").exists():
        return None

    if not worktree_is_populated():
        # A --no-checkout clone leaves HEAD at the remote tip with an empty
        # worktree. Reporting that commit would skip the checkout below.
        return None

    result = subprocess.run(
        ["git", "-C", str(TARGET), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def main() -> None:
    if shutil.which("git") is None:
        raise SystemExit("git is required but was not found on PATH")

    if TARGET.exists() and not (TARGET / ".git").exists():
        raise SystemExit(
            f"{TARGET} already exists but is not a Git checkout; move or remove it first"
        )

    if not TARGET.exists():
        TARGET.parent.mkdir(parents=True, exist_ok=True)
        run(
            "git",
            "clone",
            "--filter=blob:none",
            "--no-checkout",
            GENOPHI_REPOSITORY,
            str(TARGET),
        )
        run("git", "-C", str(TARGET), "sparse-checkout", "init", "--cone")
        run(
            "git",
            "-C",
            str(TARGET),
            "sparse-checkout",
            "set",
            "data/test_data",
        )

    if current_commit() != GENOPHI_COMMIT:
        run(
            "git",
            "-C",
            str(TARGET),
            "fetch",
            "--depth",
            "1",
            "origin",
            GENOPHI_COMMIT,
        )
        run(
            "git",
            "-C",
            str(TARGET),
            "checkout",
            "--detach",
            "FETCH_HEAD",
        )

    data_root = TARGET / "data" / "test_data"
    missing = [
        TARGET / relative
        for relative in REQUIRED_RELATIVE
        if not (TARGET / relative).exists()
    ]
    if missing:
        names = "\n".join(str(path) for path in missing)
        raise SystemExit(f"GenoPHI checkout is missing expected demo data:\n{names}")

    print(f"GenoPHI demo data ready at: {data_root}")
    print(f"Pinned commit: {GENOPHI_COMMIT}")


if __name__ == "__main__":
    main()

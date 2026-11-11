"""Bootstrap checkout tests; no network request is made in tests."""
import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/bootstrap_genophi_demo.py"
SPEC = importlib.util.spec_from_file_location("bootstrap_genophi_demo", SCRIPT)
bootstrap = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bootstrap)


def make_checkout(root, *, populated):
    """Build a fake checkout whose HEAD already equals the pinned commit."""
    target = root / "GenoPHI"
    (target / ".git").mkdir(parents=True)

    if populated:
        data_root = target / "data" / "test_data"
        (data_root / "strain_AAs").mkdir(parents=True)
        (data_root / "phage_AAs").mkdir(parents=True)
        (data_root / "ecoli_test_interaction_matrix.csv").write_text(
            "strain,phage,interaction\n"
        )

    return target


def test_empty_worktree_is_not_reported_as_the_pinned_commit(monkeypatch, tmp_path):
    """A --no-checkout clone leaves HEAD at the tip with no files on disk.

    When the pinned commit is also the remote tip, returning that commit makes
    main() skip the checkout and the demo data is never materialised.
    """
    target = make_checkout(tmp_path, populated=False)
    monkeypatch.setattr(bootstrap, "TARGET", target)

    def fail(*args):
        raise AssertionError("no subprocess should run while inspecting the checkout")

    monkeypatch.setattr(bootstrap, "run", fail)

    assert bootstrap.worktree_is_populated() is False
    assert bootstrap.current_commit() is None


def test_populated_worktree_reports_its_commit(monkeypatch, tmp_path):
    target = make_checkout(tmp_path, populated=True)
    monkeypatch.setattr(bootstrap, "TARGET", target)
    monkeypatch.setattr(
        bootstrap.subprocess,
        "run",
        lambda *args, **kwargs: type(
            "Completed", (), {"stdout": f"{bootstrap.GENOPHI_COMMIT}\n"}
        )(),
    )

    assert bootstrap.worktree_is_populated() is True
    assert bootstrap.current_commit() == bootstrap.GENOPHI_COMMIT


def test_required_paths_cover_every_demo_input():
    relative = {str(path) for path in bootstrap.REQUIRED_RELATIVE}

    assert relative == {
        "data/test_data/strain_AAs",
        "data/test_data/phage_AAs",
        "data/test_data/ecoli_test_interaction_matrix.csv",
    }

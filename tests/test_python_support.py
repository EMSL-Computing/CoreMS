"""Packaging contract for the supported Python range and TOML library."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_package_requires_python_311_and_tomlkit():
    """CoreMS installs on Python 3.11+ and declares tomlkit, not the toml package."""
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert 'requires-python = ">=3.11"' in text
    assert "tomlkit>=" in text
    assert 'target-version = "py311"' in text
    assert "Programming Language :: Python :: 3.10" not in text
    assert "Programming Language :: Python :: 3.14" in text
    for line in text.splitlines():
        stripped = line.strip()
        assert not stripped.startswith('"toml>='), stripped


def test_ci_runs_source_tests_on_python_311_and_314():
    """Regular CI tests run on Python 3.11 and 3.14 as sibling jobs.

    Notebook tests stay a single job. The two source jobs are in the same
    stage, so they start together.
    """
    text = (ROOT / ".gitlab-ci.yml").read_text(encoding="utf-8")
    assert "test-source-py311:" in text
    assert "image: python:3.11-slim" in text
    assert "test-source-py314:" in text
    assert "image: python:3.14-slim" in text
    assert "just ci-test-source" in text
    assert "test-notebooks:" in text
    assert "just ci-test-notebooks" in text

#!/usr/bin/env python3
"""
Test script to validate all example notebooks can execute without errors.
"""
import subprocess
import sys
from pathlib import Path
import argparse

# Patterns in stderr that indicate an external service is unavailable.
# Failures matching these patterns are treated as warnings (skipped) rather
# than hard failures so that transient infrastructure outages do not break CI.
EXTERNAL_SERVICE_ERROR_PATTERNS = [
    "HTTPError",
    "ConnectionError",
    "requests.exceptions",
    "503 Server Error",
    "502 Bad Gateway",
    "504 Gateway",
    "Service Temporarily Unavailable",
]


def _is_external_service_failure(stderr: str) -> bool:
    """Return True if stderr indicates an unavailable external service."""
    return any(pattern in stderr for pattern in EXTERNAL_SERVICE_ERROR_PATTERNS)


# Return values: True = pass, False = fail, None = skipped (external service)
def test_notebook(notebook_path):
    """Test a single notebook by converting it."""
    print(f"\n{'='*60}")
    print(f"Testing: {notebook_path.name}")
    print(f"{'='*60}")
    
    try:
        subprocess.run(
            [
                sys.executable,
                "-m",
                "nbconvert",
                "--to",
                "notebook",
                "--execute",
                "--ExecutePreprocessor.timeout=300",
                "--output",
                f"/tmp/{notebook_path.stem}_test.ipynb",
                str(notebook_path),
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        print(f"✓ {notebook_path.name} passed")
        return True
    except subprocess.CalledProcessError as e:
        if _is_external_service_failure(e.stderr):
            print(f"⚠ {notebook_path.name} skipped (external service unavailable)")
            print(f"STDERR:\n{e.stderr[-2000:]}")
            return None
        print(f"✗ {notebook_path.name} failed")
        print(f"STDOUT:\n{e.stdout}")
        print(f"STDERR:\n{e.stderr}")
        return False


def discover_notebooks(notebooks_dir):
    """Discover all notebooks in the notebooks directory (non-recursive)."""
    all_notebooks = []
    for notebook_path in sorted(notebooks_dir.glob("*.ipynb")):
        # Exclude checkpoint files
        if ".ipynb_checkpoints" not in str(notebook_path):
            all_notebooks.append(notebook_path)
    return all_notebooks


def main():
    """Run tests on all notebooks or a selected notebook."""
    parser = argparse.ArgumentParser(description="Execute example notebooks with nbconvert")
    parser.add_argument(
        "--notebook",
        "-n",
        help="Notebook filename or path to run (e.g., LCMS_Tutorial.ipynb)",
    )
    args = parser.parse_args()

    notebooks_dir = Path(__file__).parent / "notebooks"
    
    if not notebooks_dir.exists():
        print(f"Error: notebooks directory not found at {notebooks_dir}")
        sys.exit(1)
    
    if args.notebook:
        notebook_arg = Path(args.notebook)
        candidate = notebook_arg if notebook_arg.is_absolute() else notebooks_dir / notebook_arg
        candidate = candidate.resolve()

        if not candidate.exists():
            print(f"Error: notebook not found: {args.notebook}")
            sys.exit(1)
        notebooks = [candidate]
    else:
        # Discover all notebooks automatically
        notebooks = discover_notebooks(notebooks_dir)
    
    if not notebooks:
        print("No notebooks found to test")
        sys.exit(1)
    
    print(f"Found {len(notebooks)} notebook(s) to test")
    
    results = {}
    for notebook_path in notebooks:
        results[notebook_path.name] = test_notebook(notebook_path)
    
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    
    passed = sum(1 for v in results.values() if v is True)
    skipped = sum(1 for v in results.values() if v is None)
    failed = sum(1 for v in results.values() if v is False)
    total = len(results)
    
    for notebook, result in results.items():
        if result is True:
            status = "✓ PASS"
        elif result is None:
            status = "⚠ SKIP"
        else:
            status = "✗ FAIL"
        print(f"{status}: {notebook}")
    
    print(f"\n{passed}/{total} notebooks passed, {skipped} skipped (external service), {failed} failed")
    
    if failed > 0:
        sys.exit(1)
    
    if skipped > 0:
        print("\nSome notebooks were skipped due to unavailable external services.")
    else:
        print("\nAll tests passed!")

if __name__ == "__main__":
    main()

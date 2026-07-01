#!/usr/bin/env python3
"""
Test script to validate all example notebooks can execute without errors.
"""
import subprocess
import sys
from pathlib import Path


def test_notebook(notebook_path):
    """Test a single notebook by converting it."""
    print(f"\n{'='*60}")
    print(f"Testing: {notebook_path.name}")
    print(f"{'='*60}")
    
    try:
        result = subprocess.run(
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
    """Run tests on all notebooks."""
    notebooks_dir = Path(__file__).parent / "notebooks"
    
    if not notebooks_dir.exists():
        print(f"Error: notebooks directory not found at {notebooks_dir}")
        sys.exit(1)
    
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
    
    passed = sum(1 for v in results.values() if v)
    total = len(results)
    
    for notebook, result in results.items():
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {notebook}")
    
    print(f"\n{passed}/{total} notebooks passed")
    
    if passed < total:
        sys.exit(1)
    
    print("\nAll tests passed!")

if __name__ == "__main__":
    main()

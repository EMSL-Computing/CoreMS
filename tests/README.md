# Running CoreMS Tests

This directory contains the test suite run by `pytest`.

## Prerequisites

- Use the project virtual environment.
- Install test dependencies from the project root:

```bash
python -m pip install -r requirements.txt
python -m pip install pytest pytest-cov
```

## Common Commands

Run all tests from the project root:

```bash
pytest --cache-clear
```

Run a single test file:

```bash
pytest tests/test_wf_lipidomics.py
```

Run one test function:

```bash
pytest tests/test_wf_lipidomics.py -k lipidomics_workflow
```

## Lipidomics SQLite Test Resource

The lipidomics workflow test uses a local sqlite library file:

Download it from the project root with:

```bash
make download-lipidomics-db
```

You can override the output path when needed:

```bash
make download-lipidomics-db LIPIDOMICS_SQLITE_PATH=/path/to/202412_lipid_ref.sqlite
```

- Default local path used by tests:

```text
tests/tests_data/lcms/202412_lipid_ref.sqlite
```

- Download URL:

```text
https://nmdcdemo.emsl.pnnl.gov/minio/lipidomics/parameter_files/202412_lipid_ref.sqlite
```

By default, if the file is missing, the test fixture will try to auto-download it.

You can also provide a custom local path:

```bash
export COREMS_LIPIDOMICS_SQLITE_PATH=/path/to/202412_lipid_ref.sqlite
```

## How To Opt Out Of Lipidomics DB Tests

Skip tests that require the lipidomics sqlite database:

```bash
pytest --skip-lipidomics-db
```

You can combine this with a target file:

```bash
pytest tests/test_wf_lipidomics.py --skip-lipidomics-db
```

If you prefer not to auto-download the sqlite file during local runs, disable it:

```bash
export COREMS_LIPIDOMICS_AUTO_DOWNLOAD=0
```

With auto-download disabled, lipidomics-db tests will be skipped when the sqlite file is not present.

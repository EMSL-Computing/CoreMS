"""URL rewrite onto the psycopg 3 dialect."""

from sqlalchemy.engine.url import make_url

from corems.molecular_id.factory.postgres_url import psycopg3_url


def test_bare_postgresql_url_uses_psycopg_and_keeps_parts():
    result = psycopg3_url(
        "postgresql://user:pass@host:5432/db?sslmode=require"
    )
    parsed = make_url(result)
    assert parsed.drivername == "postgresql+psycopg"
    assert parsed.username == "user"
    assert parsed.password == "pass"
    assert parsed.host == "host"
    assert parsed.port == 5432
    assert parsed.database == "db"
    assert parsed.query["sslmode"] == "require"


def test_postgres_scheme_uses_psycopg():
    parsed = make_url(psycopg3_url("postgres://user:pass@host:5432/db"))
    assert parsed.drivername == "postgresql+psycopg"
    assert parsed.database == "db"


def test_psycopg2_url_uses_psycopg_and_keeps_encoded_password():
    parsed = make_url(
        psycopg3_url("postgresql+psycopg2://user:p%40ss@host:5432/db")
    )
    assert parsed.drivername == "postgresql+psycopg"
    assert parsed.password == "p@ss"
    assert parsed.host == "host"
    assert parsed.database == "db"


def test_psycopg3_url_is_unchanged():
    url = "postgresql+psycopg://user:pass@host/db"
    assert psycopg3_url(url) == url


def test_other_drivers_and_sentinels_are_unchanged():
    assert (
        psycopg3_url("postgresql+pg8000://user:pass@host/db")
        == "postgresql+pg8000://user:pass@host/db"
    )
    assert psycopg3_url("sqlite:///db/molformula.db") == "sqlite:///db/molformula.db"
    assert psycopg3_url(None) is None
    assert psycopg3_url("") == ""
    assert psycopg3_url("None") == "None"
    assert psycopg3_url("False") == "False"

from unittest.mock import MagicMock, patch

from corems.molecular_id.factory.EI_SQL import EI_LowRes_SQLite
from corems.molecular_id.factory.MolecularLookupTable import insert_database_worker
from corems.molecular_id.factory.molecularSQL import MolForm_SQL


def test_molform_init_engine_rewrites_psycopg2_url():
    db = MolForm_SQL.__new__(MolForm_SQL)
    with patch("corems.molecular_id.factory.molecularSQL.create_engine") as create:
        create.return_value = MagicMock()
        db.init_engine("postgresql+psycopg2://user:pass@host:5432/db")
    parsed = make_url(create.call_args.args[0])
    assert parsed.drivername == "postgresql+psycopg"
    assert parsed.password == "pass"
    assert create.call_args.kwargs["isolation_level"] == "AUTOCOMMIT"


def test_molform_init_engine_leaves_sqlite_options():
    db = MolForm_SQL.__new__(MolForm_SQL)
    with patch("corems.molecular_id.factory.molecularSQL.create_engine") as create:
        create.return_value = MagicMock()
        db.init_engine("sqlite:///db/molformula.db")
    assert create.call_args.args[0] == "sqlite:///db/molformula.db"
    assert "isolation_level" not in create.call_args.kwargs


def test_initiate_database_rewrites_bare_postgresql_url():
    db = MolForm_SQL.__new__(MolForm_SQL)
    with patch("corems.molecular_id.factory.molecularSQL.create_engine") as create:
        create.return_value.connect.return_value = MagicMock()
        db.initiate_database(
            "postgresql://user:pass@host:5432/postgres", "molformula"
        )
    assert make_url(create.call_args.args[0]).drivername == "postgresql+psycopg"


def test_ei_init_engine_rewrites_psycopg2_url():
    db = EI_LowRes_SQLite.__new__(EI_LowRes_SQLite)
    with patch("corems.molecular_id.factory.EI_SQL.create_engine") as create:
        create.return_value = MagicMock()
        db.init_engine("postgresql+psycopg2://user:pass@host:5432/lowres")
    assert make_url(create.call_args.args[0]).drivername == "postgresql+psycopg"
    assert create.call_args.kwargs["poolclass"].__name__ == "QueuePool"


def test_insert_worker_rewrites_postgres_scheme():
    with patch(
        "corems.molecular_id.factory.MolecularLookupTable.create_engine"
    ) as create, patch(
        "corems.molecular_id.factory.MolecularLookupTable.sessionmaker"
    ) as sessionmaker:
        create.return_value = MagicMock()
        sessionmaker.return_value.return_value = MagicMock()
        insert_database_worker(([], "postgres://user:pass@host:5432/db"))
    assert make_url(create.call_args.args[0]).drivername == "postgresql+psycopg"
    assert create.call_args.kwargs["isolation_level"] == "AUTOCOMMIT"

from pathlib import Path

from corems.encapsulation.factory.processingSetting import (
    MolecularFormulaSearchSettings,
)


def test_project_requires_psycopg3_sqlalchemy_21_and_python_311():
    text = Path("pyproject.toml").read_text()
    assert 'requires-python = ">=3.11"' in text
    assert "Programming Language :: Python :: 3.10" not in text
    assert "psycopg2" not in text
    assert 'psycopg[binary]>=3.2' in text
    assert "SQLAlchemy>=2.1" in text
    assert 'target-version = "py311"' in text


def test_formula_search_default_url_uses_psycopg():
    assert (
        MolecularFormulaSearchSettings.url_database
        == "postgresql+psycopg://coremsappdb:coremsapppnnl@localhost:5432/coremsapp"
    )


def test_psycopg_imports():
    import psycopg

    assert psycopg.__version__

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from sqlalchemy.engine.url import make_url

from corems.encapsulation.constant import Labels
from corems.molecular_id.factory.EI_SQL import EI_LowRes_SQLite
from corems.molecular_id.factory.MolecularLookupTable import insert_database_worker
from corems.molecular_id.factory.molecularSQL import MolForm_SQL, psycopg3_url
from corems.molecular_id.input.nistMSI import ReadNistMSI
from corems.encapsulation.factory.processingSetting  import MolecularFormulaSearchSettings

def test_nist_to_sql():

    file_location = Path.cwd() / "tests/tests_data/gcms/" / "PNNLMetV20191015.MSL"

    sqlLite_obj = ReadNistMSI(file_location).get_sqlLite_obj()

    try:
        response = sqlLite_obj.query_min_max_ri((1637.30, 1638.30)) 
        assert len(response) == 6

        response = sqlLite_obj.query_min_max_rt((17.111, 18.111))          
        assert len(response) == 137

        response = sqlLite_obj.query_min_max_ri_and_rt((1637.30, 1638.30),(17.111, 18.111)) 
        assert len(response) == 6
    finally:
        sqlLite_obj.session.close()
        sqlLite_obj.engine.dispose()

@pytest.mark.molecular_db
def test_query_sql():

    sqldb = MolForm_SQL()

    try:
        ion_type = Labels.protonated_de_ion
        classe = ['{"O": 2}']
        nominal_mz = [301]
        results = sqldb.get_dict_by_classes(classe, ion_type, nominal_mz, +1, MolecularFormulaSearchSettings())
        assert len(results.get(classe[0]).get(301)) == 3
    finally:
        sqldb.close()


@pytest.mark.parametrize(
    "url, password, query",
    [
        (
            "postgresql://user:pass@host:5432/db?sslmode=require",
            "pass",
            {"sslmode": "require"},
        ),
        ("postgres://user:pass@host:5432/db", "pass", {}),
        ("postgresql+psycopg2://user:p%40ss@host:5432/db", "p@ss", {}),
    ],
)
def test_legacy_postgres_url_uses_psycopg(url, password, query):
    parsed = make_url(psycopg3_url(url))
    assert parsed.drivername == "postgresql+psycopg"
    assert (parsed.username, parsed.password, parsed.host, parsed.port, parsed.database) == (
        "user",
        password,
        "host",
        5432,
        "db",
    )
    assert dict(parsed.query) == query


@pytest.mark.parametrize(
    "url",
    [
        "postgresql+psycopg://user:pass@host/db",
        "postgresql+pg8000://user:pass@host/db",
        "sqlite:///db/molformula.db",
        None,
        "",
        "None",
        "False",
    ],
)
def test_non_legacy_urls_are_unchanged(url):
    assert psycopg3_url(url) == url


def test_formula_engines_rewrite_postgres_and_leave_sqlite():
    mol = MolForm_SQL.__new__(MolForm_SQL)
    with patch("corems.molecular_id.factory.molecularSQL.create_engine") as create:
        create.return_value = MagicMock()
        mol.init_engine("postgresql+psycopg2://user:pass@host:5432/db")
        parsed = make_url(create.call_args.args[0])
        assert parsed.drivername == "postgresql+psycopg"
        assert parsed.password == "pass"
        assert create.call_args.kwargs["isolation_level"] == "AUTOCOMMIT"

        create.reset_mock()
        mol.init_engine("sqlite:///db/molformula.db")
        assert create.call_args.args[0] == "sqlite:///db/molformula.db"
        assert "isolation_level" not in create.call_args.kwargs

        create.reset_mock()
        create.return_value.connect.return_value = MagicMock()
        mol.initiate_database("postgresql://user:pass@host:5432/postgres", "molformula")
        assert make_url(create.call_args.args[0]).drivername == "postgresql+psycopg"

    ei = EI_LowRes_SQLite.__new__(EI_LowRes_SQLite)
    with patch("corems.molecular_id.factory.EI_SQL.create_engine") as create:
        create.return_value = MagicMock()
        ei.init_engine("postgresql+psycopg2://user:pass@host:5432/lowres")
        assert make_url(create.call_args.args[0]).drivername == "postgresql+psycopg"
        assert create.call_args.kwargs["poolclass"].__name__ == "QueuePool"

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

   
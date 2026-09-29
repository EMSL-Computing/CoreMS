"""Normalize Postgres URLs onto the psycopg 3 dialect."""

from sqlalchemy.engine.url import make_url

_PSYCOPG3_DRIVERS = {"postgresql", "postgres", "postgresql+psycopg2"}


def psycopg3_url(url):
    """Return ``url`` using the psycopg 3 dialect when it is a legacy Postgres URL.

    SQLAlchemy's default Postgres driver is psycopg 3. CoreMS does not install
    psycopg2, so a bare ``postgresql://`` or ``postgres://`` URL and a legacy
    ``postgresql+psycopg2://`` URL are rewritten to ``postgresql+psycopg://``.
    sqlite and any other explicit driver are returned unchanged.

    Parameters
    ----------
    url : str or None
        Database URL. ``None``, ``""``, ``"None"``, and ``"False"`` are
        returned as given.

    Returns
    -------
    str or None
        Rewritten URL, or ``url`` when no rewrite applies.
    """
    if url is None or url == "" or url == "None" or url == "False":
        return url
    parsed = make_url(url)
    if parsed.drivername in _PSYCOPG3_DRIVERS:
        parsed = parsed.set(drivername="postgresql+psycopg")
    return parsed.render_as_string(hide_password=False)

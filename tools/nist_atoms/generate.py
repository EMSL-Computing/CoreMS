"""Generate corems.encapsulation.nist_atoms from a vendored NIST ASCII dump.

Maintainer only. Not imported by formula search. Run via `make nist-atoms`.
Downloads the NIST dump and errors if the download fails. Writes files only
if the dump or generated tables changed.
"""

from __future__ import annotations

import re
import ssl
import sys
import urllib.error
import urllib.request
from collections import defaultdict
from dataclasses import dataclass
from datetime import date
from pathlib import Path

TOOLS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TOOLS_DIR.parents[1]
NIST_TXT = TOOLS_DIR / "AtomicWeightsAndIsotopicCompNIST.txt"
OUTPUT_PY = REPO_ROOT / "corems" / "encapsulation" / "nist_atoms.py"
CHANGES_MD = TOOLS_DIR / "CHANGES.md"

NIST_LANDING_URL = (
    "https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses"
)
NIST_DUMP_URL = (
    "https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl"
    "?ele=&ascii=ascii2&isotype=all"
)
NIST_TABLE_ID = (
    "Coursey et al. Atomic Weights and Isotopic Compositions version 4.1. "
    + NIST_LANDING_URL
)

MASS_SIGNIFICANT = 1e-7
ABUNDANCE_SIGNIFICANT = 1e-6

_UNCERT = re.compile(r"\([^)]*\)")

# CHONPS / group order for bare element symbols (formula .string stability).
ELEMENT_ORDER = [
    "C", "H", "O", "N", "P", "S",
    "F", "Cl", "Br", "I",
    "Li", "Na", "K", "Rb", "Cs",
    "He", "Ne", "Ar", "Kr", "Xe",
    "Be", "B",
    "Mg", "Al", "Si",
    "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
    "Ge", "As", "Se",
    "Sr", "Y", "Zr", "Nb", "Mo", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te",
    "Ba", "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi",
    "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Th", "Pa", "U",
]

HYDROGEN_CANONICAL = {1: "H", 2: "D", 3: "T"}


@dataclass(frozen=True)
class Nuclide:
    z: int
    nist_symbol: str
    element: str
    mass_number: int
    mass: float
    mass_literal: str
    abundance: float
    abundance_literal: str


def _field_value(line: str) -> str:
    return line.split("=", 1)[1].strip()


def parse_nist_number(raw: str) -> tuple[float | None, str | None]:
    cleaned = raw.strip().replace("&nbsp;", "").replace("#", "")
    cleaned = _UNCERT.sub("", cleaned).strip()
    if not cleaned:
        return None, None
    return float(cleaned), cleaned


def _read_header_retrieved(text: str) -> str:
    for line in text.splitlines():
        if line.startswith("# Retrieved:"):
            return line.split(":", 1)[1].strip()
    return ""


def records_body(text: str) -> str:
    """Return normalized NIST 8-line records, stripping HTML chrome and comments."""
    lowered = text.lower()
    start = lowered.find("<pre>")
    end = lowered.find("</pre>")
    if start != -1 and end != -1:
        text = text[start + len("<pre>") : end]
    text = text.replace("&nbsp;", "")
    lines: list[str] = []
    started = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        if not started:
            if stripped.startswith("Atomic Number"):
                started = True
            else:
                continue
        lines.append(line.rstrip())
    while lines and not lines[-1].strip():
        lines.pop()
    if not started:
        return ""
    return "\n".join(lines) + "\n"


def _ssl_context():
    try:
        import certifi
    except ImportError:
        return None
    return ssl.create_default_context(cafile=certifi.where())


def fetch_nist_dump() -> str:
    """Download the linearized ASCII all-isotopes dump. Errors if unavailable."""
    request = urllib.request.Request(
        NIST_DUMP_URL,
        headers={"User-Agent": "CoreMS nist-atoms generator"},
    )
    context = _ssl_context()
    try:
        with urllib.request.urlopen(request, timeout=60, context=context) as response:
            raw = response.read()
            status = getattr(response, "status", None)
            if status is not None and status >= 400:
                raise SystemExit(
                    f"Failed to download NIST dump from {NIST_DUMP_URL}: HTTP {status}"
                )
    except urllib.error.URLError as exc:
        hint = ""
        if "CERTIFICATE" in str(exc).upper():
            hint = (
                " Use the CoreMS venv Python (certifi), e.g. "
                "PYTHON=.venv/bin/python make nist-atoms."
            )
        raise SystemExit(
            f"Failed to download NIST dump from {NIST_DUMP_URL}: {exc}.{hint}"
        ) from exc
    text = raw.decode("utf-8", errors="replace")
    if "Atomic Number" not in text:
        raise SystemExit(
            f"NIST dump from {NIST_DUMP_URL} did not contain Atomic Number records"
        )
    return text


def vendored_header(retrieved: str) -> str:
    return (
        "# NIST Atomic Weights and Isotopic Compositions "
        "(linearized ASCII, all isotopes).\n"
        f"# Landing: {NIST_LANDING_URL}\n"
        f"# Dump: {NIST_DUMP_URL}\n"
        f"# Retrieved: {retrieved}\n"
        "# Compilation: Coursey et al. version 4.1 (https://physics.nist.gov/Comp)\n"
        "#\n"
        "# Do not edit records by hand. `make nist-atoms` re-downloads this file.\n"
        "# Parser ignores lines starting with #.\n"
        "\n"
    )


def parse_nist(path: Path) -> tuple[list[Nuclide], str]:
    text = path.read_text(encoding="utf-8")
    retrieved = _read_header_retrieved(text)
    records: list[dict[str, str]] = []
    current: dict[str, str] = {}
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            if current and "Atomic Number" in current:
                records.append(current)
                current = {}
            continue
        if "=" not in line:
            continue
        key, _sep, _rest = line.partition("=")
        key = key.strip()
        if key == "Atomic Number" and current:
            records.append(current)
            current = {}
        current[key] = _field_value(line)
    if current and "Atomic Number" in current:
        records.append(current)

    nuclides: list[Nuclide] = []
    for rec in records:
        z = int(rec["Atomic Number"])
        nist_symbol = rec["Atomic Symbol"].strip()
        mass_number = int(rec["Mass Number"])
        mass, mass_lit = parse_nist_number(rec.get("Relative Atomic Mass", ""))
        abundance, abund_lit = parse_nist_number(rec.get("Isotopic Composition", ""))
        if mass is None:
            continue
        if abundance is None:
            continue
        element = "H" if z == 1 else nist_symbol
        nuclides.append(
            Nuclide(
                z=z,
                nist_symbol=nist_symbol,
                element=element,
                mass_number=mass_number,
                mass=mass,
                mass_literal=mass_lit,
                abundance=abundance,
                abundance_literal=abund_lit,
            )
        )
    return nuclides, retrieved


def canonical_key(n: Nuclide, most_abundant: Nuclide) -> str:
    if n.z == 1:
        return HYDROGEN_CANONICAL[n.mass_number]
    if n.mass_number == most_abundant.mass_number:
        return n.element
    return f"{n.mass_number}{n.element}"


def nuclide_lookup_key(n: Nuclide) -> str:
    return f"{n.mass_number}{n.element}"


def build_tables(nuclides: list[Nuclide], retrieved: str) -> dict:
    by_element: dict[str, list[Nuclide]] = defaultdict(list)
    for n in nuclides:
        by_element[n.element].append(n)

    atomic_masses: dict[str, float] = {}
    mass_literals: dict[str, str] = {}
    isotopic_abundance: dict[str, float] = {}
    abund_literals: dict[str, str] = {}
    isotopes: dict[str, list] = {}
    canonical_of: dict[tuple[str, int], str] = {}

    for element, group in by_element.items():
        most = max(group, key=lambda n: n.abundance)
        rares: list[tuple[int, str]] = []
        for n in group:
            ckey = canonical_key(n, most)
            nkey = nuclide_lookup_key(n)
            canonical_of[(element, n.mass_number)] = ckey
            atomic_masses[ckey] = n.mass
            mass_literals[ckey] = n.mass_literal
            isotopic_abundance[ckey] = n.abundance
            abund_literals[ckey] = n.abundance_literal
            if nkey != ckey:
                atomic_masses[nkey] = n.mass
                mass_literals[nkey] = n.mass_literal
                isotopic_abundance[nkey] = n.abundance
                abund_literals[nkey] = n.abundance_literal
            most_canonical = canonical_key(most, most)
            if ckey != most_canonical:
                rares.append((n.mass_number, ckey))
        rares.sort()
        rare_keys = [key for _a, key in rares]
        isotopes[element] = rare_keys if rare_keys else [None]

    present = set(isotopes)
    extra = sorted(present - set(ELEMENT_ORDER))
    element_block = [e for e in ELEMENT_ORDER if e in present] + extra
    rare_tail = []
    for element in element_block:
        rare_keys = isotopes[element]
        if rare_keys and rare_keys[0] is not None:
            rare_tail.extend(rare_keys)
    rare_tail.sort(key=lambda k: (_mass_number_from_key(k), k))
    atoms_order = element_block + rare_tail

    return {
        "NIST_TABLE_ID": NIST_TABLE_ID,
        "NIST_DUMP_URL": NIST_DUMP_URL,
        "NIST_RETRIEVED": retrieved,
        "atomic_masses": atomic_masses,
        "isotopic_abundance": isotopic_abundance,
        "isotopes": isotopes,
        "atoms_order": atoms_order,
        "mass_literals": mass_literals,
        "abund_literals": abund_literals,
        "canonical_keys": [k for k in atoms_order],
    }


def _mass_number_from_key(key: str) -> int:
    digits = "".join(ch for ch in key if ch.isdigit())
    return int(digits) if digits else 0


def _format_str_dict(mapping: dict[str, str], key_order: list[str]) -> str:
    lines = ["{"]
    seen = set()
    for key in key_order:
        if key in mapping:
            lines.append(f"    {key!r}: {mapping[key]},")
            seen.add(key)
    for key in sorted(mapping):
        if key not in seen:
            lines.append(f"    {key!r}: {mapping[key]},")
    lines.append("}")
    return "\n".join(lines)


def _format_isotopes(isotopes: dict[str, list], element_order: list[str]) -> str:
    lines = ["{"]
    for element in element_order:
        if element not in isotopes:
            continue
        rares = isotopes[element]
        lines.append(f"    {element!r}: {rares!r},")
    lines.append("}")
    return "\n".join(lines)


def _format_list(values: list[str]) -> str:
    lines = ["["]
    for item in values:
        lines.append(f"    {item!r},")
    lines.append("]")
    return "\n".join(lines)


def render_module(tables: dict) -> str:
    key_order = tables["canonical_keys"]
    body = f'''"""NIST-pinned atomic masses, abundances, and isotope lists.

Do not edit by hand. Regenerate with `make nist-atoms`.

Hydrogen aliases: H/1H, D/2H.
Most-abundant nuclides are stored under both the bare symbol and the
mass-number key (C and 12C). Formula strings use canonical keys only.

``isotopes`` maps an element symbol to rare-nuclide keys, or ``[None]``
when there is no heavy isotope (same sentinel ``Atoms`` already uses).
English names are not a NIST field; see ``Atoms.element_names``.
"""

NIST_TABLE_ID = {tables["NIST_TABLE_ID"]!r}
NIST_DUMP_URL = {tables["NIST_DUMP_URL"]!r}
NIST_RETRIEVED = {tables["NIST_RETRIEVED"]!r}

atomic_masses = {_format_str_dict(tables["mass_literals"], key_order)}

isotopic_abundance = {_format_str_dict(tables["abund_literals"], key_order)}

isotopes = {_format_isotopes(tables["isotopes"], [e for e in tables["atoms_order"] if e in tables["isotopes"]])}

atoms_order = {_format_list(tables["atoms_order"])}
'''
    return body


def load_previous_tables() -> tuple[dict[str, float], dict[str, float], set[str], str]:
    if OUTPUT_PY.is_file():
        ns: dict = {}
        exec(OUTPUT_PY.read_text(encoding="utf-8"), ns)
        return (
            dict(ns["atomic_masses"]),
            dict(ns["isotopic_abundance"]),
            set(ns["atoms_order"]),
            ns.get("NIST_TABLE_ID", "previous nist_atoms.py"),
        )
    sys.path.insert(0, str(REPO_ROOT))
    from corems.encapsulation.constant import Atoms  # noqa: WPS433

    return (
        dict(Atoms.atomic_masses),
        dict(Atoms.isotopic_abundance),
        set(Atoms.atoms_order),
        "corems.encapsulation.constant.Atoms (pre-NIST module)",
    )


def write_changes(
    tables: dict,
    previous_masses: dict,
    previous_abund: dict,
    previous_canonical: set[str],
    previous_id: str,
) -> None:
    new_masses: dict[str, float] = tables["atomic_masses"]
    new_abund: dict[str, float] = tables["isotopic_abundance"]
    new_canon = set(tables["canonical_keys"])
    old_canon = set(previous_canonical)

    added = sorted(new_canon - old_canon)
    removed = sorted(old_canon - new_canon)
    shared = sorted(new_canon & old_canon)

    all_rows = []
    significant = []
    for key in shared:
        old_m = previous_masses.get(key)
        new_m = new_masses.get(key)
        old_a = previous_abund.get(key)
        new_a = new_abund.get(key)
        d_m = None if old_m is None or new_m is None else new_m - old_m
        d_a = None if old_a is None or new_a is None else new_a - old_a
        if (d_m is None or d_m == 0) and (d_a is None or d_a == 0):
            continue
        all_rows.append((key, d_m, d_a, old_m, new_m, old_a, new_a))
        sig_m = d_m is not None and abs(d_m) >= MASS_SIGNIFICANT
        sig_a = d_a is not None and abs(d_a) >= ABUNDANCE_SIGNIFICANT
        if sig_m or sig_a:
            significant.append((key, d_m, d_a, old_m, new_m, old_a, new_a))

    def fmt_delta(key, d_m, d_a, old_m, new_m, old_a, new_a) -> str:
        bits = [f"- `{key}`"]
        if d_m:
            bits.append(f"mass {old_m} → {new_m} (Δ {d_m:+.8g} u)")
            if abs(d_m) >= 0.5:
                bits.append(
                    "(most-abundant nuclide assignment likely changed)"
                )
        if d_a:
            bits.append(f"abundance {old_a} → {new_a} (Δ {d_a:+.8g})")
        return " ".join(bits)

    lines = [
        "# NIST atoms change log",
        "",
        f"Previous: {previous_id}",
        f"New: {tables['NIST_TABLE_ID']} (retrieved {tables['NIST_RETRIEVED']})",
        "",
        "## Added canonical keys",
        "",
    ]
    lines.extend(f"- `{k}`" for k in added) if added else lines.append("- None")
    lines += ["", "## Removed canonical keys", ""]
    lines.extend(f"- `{k}`" for k in removed) if removed else lines.append("- None")
    lines += ["", "## Significant (copy into release notes)", ""]
    if significant or added or removed:
        if added:
            lines.append("Added: " + ", ".join(f"`{k}`" for k in added))
        if removed:
            lines.append("Removed: " + ", ".join(f"`{k}`" for k in removed))
        lines.extend(fmt_delta(*row) for row in significant)
        if not significant and (added or removed):
            lines.append("No mass/abundance deltas above threshold; see added/removed above.")
    else:
        lines.append("None")
    lines += ["", "## All mass/abundance deltas (canonical keys)", ""]
    if all_rows:
        lines.extend(fmt_delta(*row) for row in all_rows)
    else:
        lines.append("None")
    lines.append("")
    CHANGES_MD.write_text("\n".join(lines), encoding="utf-8")


def check_committed_module() -> None:
    nuclides, retrieved = parse_nist(NIST_TXT)
    tables = build_tables(nuclides, retrieved)
    if not OUTPUT_PY.is_file():
        raise SystemExit(f"missing {OUTPUT_PY}")
    ns: dict = {}
    exec(OUTPUT_PY.read_text(encoding="utf-8"), ns)
    for name in ("atomic_masses", "isotopic_abundance", "isotopes", "atoms_order"):
        if ns[name] != tables[name]:
            raise SystemExit(f"nist_atoms.py is out of date ({name} mismatch). Run make nist-atoms.")
    if ns.get("NIST_TABLE_ID") != tables["NIST_TABLE_ID"]:
        raise SystemExit("nist_atoms.py NIST_TABLE_ID mismatch. Run make nist-atoms.")
    if ns.get("NIST_DUMP_URL") != tables["NIST_DUMP_URL"]:
        raise SystemExit("nist_atoms.py NIST_DUMP_URL mismatch. Run make nist-atoms.")


def generate() -> None:
    fetched = fetch_nist_dump()
    new_body = records_body(fetched)
    if not new_body:
        raise SystemExit("NIST dump contained no Atomic Number records after HTML strip")

    existing_text = NIST_TXT.read_text(encoding="utf-8") if NIST_TXT.is_file() else ""
    existing_body = records_body(existing_text) if existing_text else ""
    dump_changed = new_body != existing_body

    if dump_changed:
        retrieved = date.today().isoformat()
        NIST_TXT.write_text(vendored_header(retrieved) + new_body, encoding="utf-8")
        print(f"Updated {NIST_TXT.relative_to(REPO_ROOT)}")
    else:
        retrieved = _read_header_retrieved(existing_text)
        print("NIST dump unchanged")

    # Parse the fetched or existing NIST into nuclides
    nuclides, _header_retrieved = parse_nist(NIST_TXT)
    tables = build_tables(nuclides, retrieved or _header_retrieved)
    new_module = render_module(tables)
    existing_module = OUTPUT_PY.read_text(encoding="utf-8") if OUTPUT_PY.is_file() else ""
    if not dump_changed and existing_module == new_module:
        print("nist_atoms.py already up to date; no files written")
        return

    previous_masses, previous_abund, previous_canonical, previous_id = (
        load_previous_tables()
    )
    OUTPUT_PY.write_text(new_module, encoding="utf-8")
    write_changes(
        tables, previous_masses, previous_abund, previous_canonical, previous_id
    )
    print(f"Wrote {OUTPUT_PY.relative_to(REPO_ROOT)}")
    print(f"Wrote {CHANGES_MD.relative_to(REPO_ROOT)}")


def main() -> int:
    generate()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

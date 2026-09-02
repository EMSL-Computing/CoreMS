"""Generate corems.encapsulation.nist_atoms from a vendored NIST ASCII dump.

Maintainer only. Run via `make nist-atoms`. Downloads the NIST dump and errors
if that fails. Writes files only if the dump or generated tables changed.
"""

from __future__ import annotations

import re
import ssl
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
    element: str
    mass_number: int
    mass: float
    mass_literal: str
    abundance: float
    abundance_literal: str


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
    """Normalized NIST 8-line records, HTML chrome and comments stripped."""
    lowered = text.lower()
    start, end = lowered.find("<pre>"), lowered.find("</pre>")
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
    return ("\n".join(lines) + "\n") if started else ""


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
    try:
        with urllib.request.urlopen(
            request, timeout=60, context=_ssl_context()
        ) as response:
            raw = response.read()
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


def _nuclide_from_record(rec: dict[str, str]) -> Nuclide | None:
    mass, mass_lit = parse_nist_number(rec.get("Relative Atomic Mass", ""))
    abundance, abund_lit = parse_nist_number(rec.get("Isotopic Composition", ""))
    if mass is None or abundance is None:
        return None
    z = int(rec["Atomic Number"])
    nist_symbol = rec["Atomic Symbol"].strip()
    return Nuclide(
        z=z,
        element="H" if z == 1 else nist_symbol,
        mass_number=int(rec["Mass Number"]),
        mass=mass,
        mass_literal=mass_lit,
        abundance=abundance,
        abundance_literal=abund_lit,
    )


def parse_nist(path: Path) -> tuple[list[Nuclide], str]:
    text = path.read_text(encoding="utf-8")
    nuclides: list[Nuclide] = []
    rec: dict[str, str] = {}
    for line in records_body(text).splitlines():
        if not line.strip():
            n = _nuclide_from_record(rec)
            if n:
                nuclides.append(n)
            rec = {}
            continue
        if "=" in line:
            key, _, value = line.partition("=")
            rec[key.strip()] = value.strip()
    n = _nuclide_from_record(rec)
    if n:
        nuclides.append(n)
    return nuclides, _read_header_retrieved(text)


def canonical_key(n: Nuclide, most_abundant: Nuclide) -> str:
    if n.z == 1:
        return HYDROGEN_CANONICAL[n.mass_number]
    if n.mass_number == most_abundant.mass_number:
        return n.element
    return f"{n.mass_number}{n.element}"


def nuclide_lookup_key(n: Nuclide) -> str:
    return f"{n.mass_number}{n.element}"


def _store(n: Nuclide, key: str, masses: dict, mass_lits: dict, abunds: dict, abund_lits: dict) -> None:
    masses[key] = n.mass
    mass_lits[key] = n.mass_literal
    abunds[key] = n.abundance
    abund_lits[key] = n.abundance_literal


def build_tables(nuclides: list[Nuclide], retrieved: str) -> dict:
    by_element: dict[str, list[Nuclide]] = defaultdict(list)
    for n in nuclides:
        by_element[n.element].append(n)

    atomic_masses: dict[str, float] = {}
    mass_literals: dict[str, str] = {}
    isotopic_abundance: dict[str, float] = {}
    abund_literals: dict[str, str] = {}
    isotopes: dict[str, list] = {}

    for element, group in by_element.items():
        most = max(group, key=lambda n: n.abundance)
        most_canonical = canonical_key(most, most)
        rares: list[tuple[int, str]] = []
        for n in group:
            ckey = canonical_key(n, most)
            _store(n, ckey, atomic_masses, mass_literals, isotopic_abundance, abund_literals)
            nkey = nuclide_lookup_key(n)
            if nkey != ckey:
                _store(n, nkey, atomic_masses, mass_literals, isotopic_abundance, abund_literals)
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
        "canonical_keys": list(atoms_order),
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
        if element in isotopes:
            lines.append(f"    {element!r}: {isotopes[element]!r},")
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
    iso_order = [e for e in tables["atoms_order"] if e in tables["isotopes"]]
    return f'''"""NIST-pinned atomic masses, abundances, and isotope lists.

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

isotopes = {_format_isotopes(tables["isotopes"], iso_order)}

atoms_order = {_format_list(tables["atoms_order"])}
'''


def _exec_module(path: Path) -> dict:
    ns: dict = {}
    exec(path.read_text(encoding="utf-8"), ns)
    return ns


def load_previous_tables() -> tuple[dict[str, float], dict[str, float], set[str], str]:
    ns = _exec_module(OUTPUT_PY)
    return (
        dict(ns["atomic_masses"]),
        dict(ns["isotopic_abundance"]),
        set(ns["atoms_order"]),
        ns.get("NIST_TABLE_ID", "previous nist_atoms.py"),
    )


def write_changes(
    tables: dict,
    previous_masses: dict,
    previous_abund: dict,
    previous_canonical: set[str],
    previous_id: str,
) -> None:
    new_masses = tables["atomic_masses"]
    new_abund = tables["isotopic_abundance"]
    new_canon = set(tables["canonical_keys"])
    old_canon = set(previous_canonical)
    added = sorted(new_canon - old_canon)
    removed = sorted(old_canon - new_canon)

    all_rows = []
    significant = []
    for key in sorted(new_canon & old_canon):
        old_m, new_m = previous_masses.get(key), new_masses.get(key)
        old_a, new_a = previous_abund.get(key), new_abund.get(key)
        d_m = None if old_m is None or new_m is None else new_m - old_m
        d_a = None if old_a is None or new_a is None else new_a - old_a
        if not d_m and not d_a:
            continue
        row = (key, d_m, d_a, old_m, new_m, old_a, new_a)
        all_rows.append(row)
        if (d_m and abs(d_m) >= MASS_SIGNIFICANT) or (
            d_a and abs(d_a) >= ABUNDANCE_SIGNIFICANT
        ):
            significant.append(row)

    def fmt_delta(key, d_m, d_a, old_m, new_m, old_a, new_a) -> str:
        bits = [f"- `{key}`"]
        if d_m:
            bits.append(f"mass {old_m} → {new_m} (Δ {d_m:+.8g} u)")
            if abs(d_m) >= 0.5:
                bits.append("(most-abundant nuclide assignment likely changed)")
        if d_a:
            bits.append(f"abundance {old_a} → {new_a} (Δ {d_a:+.8g})")
        return " ".join(bits)

    def bullets(keys: list[str]) -> list[str]:
        return [f"- `{k}`" for k in keys] if keys else ["- None"]

    lines = [
        "# NIST atoms change log",
        "",
        f"Previous: {previous_id}",
        f"New: {tables['NIST_TABLE_ID']} (retrieved {tables['NIST_RETRIEVED']})",
        "",
        "## Added canonical keys",
        "",
        *bullets(added),
        "",
        "## Removed canonical keys",
        "",
        *bullets(removed),
        "",
        "## Significant (copy into release notes)",
        "",
    ]
    if significant or added or removed:
        if added:
            lines.append("Added: " + ", ".join(f"`{k}`" for k in added))
        if removed:
            lines.append("Removed: " + ", ".join(f"`{k}`" for k in removed))
        lines.extend(fmt_delta(*row) for row in significant)
        if not significant and (added or removed):
            lines.append(
                "No mass/abundance deltas above threshold; see added/removed above."
            )
    else:
        lines.append("None")
    lines += ["", "## All mass/abundance deltas (canonical keys)", ""]
    lines.extend(fmt_delta(*row) for row in all_rows) if all_rows else lines.append("None")
    lines.append("")
    CHANGES_MD.write_text("\n".join(lines), encoding="utf-8")


def check_committed_module() -> None:
    nuclides, retrieved = parse_nist(NIST_TXT)
    tables = build_tables(nuclides, retrieved)
    if not OUTPUT_PY.is_file():
        raise SystemExit(f"missing {OUTPUT_PY}")
    ns = _exec_module(OUTPUT_PY)
    for name in ("atomic_masses", "isotopic_abundance", "isotopes", "atoms_order"):
        if ns[name] != tables[name]:
            raise SystemExit(
                f"nist_atoms.py is out of date ({name} mismatch). Run make nist-atoms."
            )
    if ns.get("NIST_TABLE_ID") != tables["NIST_TABLE_ID"]:
        raise SystemExit("nist_atoms.py NIST_TABLE_ID mismatch. Run make nist-atoms.")
    if ns.get("NIST_DUMP_URL") != tables["NIST_DUMP_URL"]:
        raise SystemExit("nist_atoms.py NIST_DUMP_URL mismatch. Run make nist-atoms.")


def generate() -> None:
    new_body = records_body(fetch_nist_dump())
    if not new_body:
        raise SystemExit("NIST dump contained no Atomic Number records after HTML strip")

    existing_text = NIST_TXT.read_text(encoding="utf-8") if NIST_TXT.is_file() else ""
    dump_changed = new_body != records_body(existing_text)

    if dump_changed:
        retrieved = date.today().isoformat()
        NIST_TXT.write_text(vendored_header(retrieved) + new_body, encoding="utf-8")
        print(f"Updated {NIST_TXT.relative_to(REPO_ROOT)}")
    else:
        retrieved = _read_header_retrieved(existing_text)
        print("NIST dump unchanged")

    nuclides, header_retrieved = parse_nist(NIST_TXT)
    tables = build_tables(nuclides, retrieved or header_retrieved)
    new_module = render_module(tables)
    existing_module = OUTPUT_PY.read_text(encoding="utf-8") if OUTPUT_PY.is_file() else ""
    if not dump_changed and existing_module == new_module:
        print("nist_atoms.py already up to date; no files written")
        return

    previous = load_previous_tables()
    OUTPUT_PY.write_text(new_module, encoding="utf-8")
    write_changes(tables, *previous)
    print(f"Wrote {OUTPUT_PY.relative_to(REPO_ROOT)}")
    print(f"Wrote {CHANGES_MD.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    generate()

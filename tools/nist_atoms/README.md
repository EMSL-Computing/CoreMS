# NIST atom tables (maintainer)

CoreMS ships a **parsed, pinned** NIST snapshot. Formula search only looks up
`Atoms.atomic_masses` / abundances — it never parses NIST at runtime.

## Source

Human page:

https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses

Download settings: All Elements, Linearized ASCII, All isotopes.

Dump URL:

https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=all

`make nist-atoms` downloads that dump (errors if the download fails), strips HTML
chrome, and compares records to the vendored file. If the dump is unchanged and
`nist_atoms.py` already matches, it writes nothing.

## Commands

```bash
make nist-atoms   # download NIST dump; rewrite tables only if needed
```

Do not hook this into `make patch|minor|major`. Release prep runs `make nist-atoms`
**before** `make lint`.

If the pin changed, paste the **Breaking** and **Significant** sections of
`CHANGES.md` into release notes (`RELEASE.md`). Significant lists **Deleted
from the library** (no mass lookup) separately from **lookup aliases** that
remain on `Atoms.atomic_masses` but are not formula-string keys. Do not
describe alias keys (`114Cd`, `12C`, `40Ca`, …) as removed.

## Policy

- Include every nuclide NIST lists with an isotopic composition; skip empty composition (no T, ¹⁴C allowlist)
- Canonical formula keys: most abundant = `C`; rares = `13C`; hydrogen = `H`/`D`
- Lookup aliases on masses/abundances: `12C`, `1H`, `2H`
- English names are **not** NIST data; they live on `Atoms.element_names` in `constant.py`
- IsoSpec is the isotopologue enumerator only; it is still passed CoreMS/NIST numbers

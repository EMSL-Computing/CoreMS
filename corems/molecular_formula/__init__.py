"""Molecular formula objects and related calculations.

This package defines molecular formula types used throughout CoreMS
(for example after molecular formula search) and calculations such as
exact mass, double-bond equivalent, and assignment confidence.

**Confidence score**

When a formula is attached to a mass spectral peak (and isotopologues have
been evaluated during search), CoreMS provides a composite **confidence
score** that ranks candidates. Public properties on
``MolecularFormula``
(:class:`~corems.molecular_formula.factory.MolecularFormulaFactory.MolecularFormula`)
are:

- ``mz_error_score`` — mass accuracy of the monoisotopic peak (Gaussian
  score from the ppm error).
- ``average_mz_error_score`` — mean mass-accuracy score over the mono peak
  and expected isotopologue peaks (missing expected isotopologues score 0).
- ``isotopologue_similarity`` — agreement between expected and observed
  isotopologue abundance patterns (Manhattan-based similarity).
- ``confidence_score`` — weighted sum of the formula-level mass-error term
  and the isotopologue similarity term.

**Composite score**

```
CS = (w_err * m_err) + (w_iso * m_iso)
```

where:

- ``m_err`` is ``average_mz_error_score``
- ``m_iso`` is ``isotopologue_similarity``
- ``w_err`` is ``molecular_search_settings.mz_error_score_weight``
  (default **0.6**)
- ``w_iso`` is ``molecular_search_settings.isotopologue_score_weight``
  (default **0.4**)

**Mass-error term (per peak)**

For a single peak, ``mz_error_score`` is a Gaussian score of the assignment
error ``delta`` (ppm) about a mean of zero (calibrated spectrum assumed):

```
m = exp( -(delta ** 2) / (2 * (sigma ** 2)) )
```

The width ``sigma`` is the peak's ``predicted_std`` when set (for example from
resolving-power / mass-error prediction). If unset, CoreMS uses a fallback
of **1.66** ppm. Formula-level ``average_mz_error_score`` averages this score
over the monoisotopic assignment and each expected isotopologue peak.

**Isotopologue term**

``isotopologue_similarity`` compares sum-normalized theoretical and
observed abundances for the mono peak plus expected isotopologues. Missing
expected partners are treated as near-zero abundance. If no isotopologues
are expected (or none were computed), the similarity is **0**.

**Interpretation**

Scores are typically in about [0, 1], with higher values preferred.
When ``score_method`` is ``"prob_score"``, ranking uses this composite
confidence score. Isotopologue detection during search
(``find_isotopologues=True``) is required for a meaningful
``isotopologue_similarity`` contribution.

**Literature**

The confidence score formulation and its use for ranking formula
candidates are described in Dewey, Corilo, Kew, and Boiteau,
*Analytical Chemistry* **2025**, *97*, 13031–13039
(https://doi.org/10.1021/acs.analchem.4c06826).

**See also**

- ``MolecularFormula.confidence_score``
- ``MolecularFormula.mz_error_score``
- ``MolecularFormula.average_mz_error_score``
- ``MolecularFormula.isotopologue_similarity``
- Settings: ``mz_error_score_weight``, ``isotopologue_score_weight``
"""

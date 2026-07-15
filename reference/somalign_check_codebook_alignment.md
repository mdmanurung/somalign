# Check query-reference codebook alignment before cell-level computation

Compares a query SOM codebook (already in reference-scaled space)
against the reference codebook across three diagnostic dimensions. The
check is O(nodes\\^2 \times\\ p) — milliseconds on a 900-node SOM — and
is designed to surface distribution mismatches before any per-cell work
begins.

## Usage

``` r
somalign_check_codebook_alignment(
  query_codebook,
  reference,
  query_masses = NULL,
  epsilon = 0.1,
  stop_if_critical = TRUE
)
```

## Arguments

- query_codebook:

  Numeric matrix of query SOM codebook vectors in reference-scaled
  coordinate space (nodes \\\times\\ features). Column names must
  include all features in `reference$features`.

- reference:

  A `somalign_reference` object.

- query_masses:

  Optional numeric vector of query node masses (length
  `nrow(query_codebook)`). Used for the mass-weighted centroid check.
  When `NULL`, uniform weights are assumed.

- epsilon:

  The OT regularisation parameter that will be passed to
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).
  Used to contextualise the cost matrix coverage check. Default `0.1`.

- stop_if_critical:

  If `TRUE` (default), throw an error when any feature has zero range
  overlap. Set to `FALSE` to emit a warning instead and return the
  diagnostics.

## Value

A `somalign_codebook_check` list (returned invisibly) with:

- `per_feature`:

  Data frame: one row per feature with `ref_min`, `ref_max`,
  `query_min`, `query_max`, `overlap_fraction`, `centroid_drift`,
  `centroid_drift_sd`, and `flag` (`"ok"` / `"warning"` / `"critical"`).

- `cost_summary`:

  Named numeric vector: `median_cost` (raw median pairwise squared
  distance), `p95_cost`, `cost_scale` (normalisation factor used by
  `somalign_fit`), and `fraction_near_eps` (fraction of pairs within
  \\3\varepsilon\\ of the normalised cost).

- `n_critical_features`:

  Number of features with zero overlap.

- `n_warning_features`:

  Number of features flagged as warning.

- `verdict`:

  `"pass"`, `"warning"`, or `"critical"`.

## Details

**Range overlap (per feature):** Does the query codebook's value range
for each marker intersect the reference codebook's range? Zero overlap
means every query node sits entirely outside the reference for that
marker: a critical failure. Less than 50\\

**Mass-weighted centroid drift (per feature):** The mass-weighted mean
of the query codebook minus that of the reference, expressed in units of
the reference codebook standard deviation. Drift \> 3 SDs flags a global
batch shift that the OT plan may not be able to absorb.

**Transport coverage (cost matrix preview):** Fraction of
query-reference codebook pairs whose normalised squared distance falls
within \\3\varepsilon\\. Pairs outside this band contribute negligible
weight to the Sinkhorn kernel. If fewer than 1\\ \\3\varepsilon\\, the
transport plan will be near-singular and most query mass will be
destroyed.

## See also

[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md),
[`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md)

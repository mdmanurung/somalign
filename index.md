# somalign

`somalign` aligns query self-organising maps to fixed `kohonen`
reference maps using codebook-level unbalanced entropic optimal
transport.

The conservative primary result is direct projection into the reference
SOM. Transport-corrected projections are returned as auxiliary columns
for visualisation, annotation, and triage.

## Installation

``` r

# install.packages("pak")
pak::pkg_install("mdmanurung/somalign")
```

## Quick start

``` r

library(kohonen)
library(somalign)

# Train reference SOM on old/reference data
reference <- somalign_train_reference(old_matrix, labels = old_labels)

# Build query object from new data
query <- somalign_query(new_matrix, reference)

# Align and extract per-sample results
fit     <- somalign_fit(query, reference)
results <- somalign_results(fit)
```

Key output columns in `results`:

| Column | Description |
|----|----|
| `old_som_unit` | Direct reference-node assignment (primary result) |
| `old_som_distance` | Distance to assigned reference node |
| `outside_reference_distance` | `TRUE` if distance exceeds reference quantile threshold |
| `final_status` | `"inside_reference"` / `"outside_reference"` / `"unknown_reference_distance"` |
| `old_som_label` | Majority label of the assigned reference node |
| `old_som_label_confidence` | Fraction of that node’s training mass with the majority label |
| `corrected_som_unit` | OT-corrected reference-node assignment (auxiliary) |
| `correction_norm` | Magnitude of the OT correction shift |
| `transferred_label` | Label transferred via OT correspondence (auxiliary) |
| `transferred_label_confidence` | Confidence of the transferred label |
| `transferred_label_accepted` | `TRUE` if match fraction and confidence thresholds are met |

## Learn more

- [`vignette("somalign", package = "somalign")`](https://mdmanurung.github.io/somalign/articles/somalign.md):
  quick start from reference training to result inspection.
- [`vignette("pretrained-old-and-new-soms", package = "somalign")`](https://mdmanurung.github.io/somalign/articles/pretrained-old-and-new-soms.md):
  existing SOMs, saved codebooks, diagnostics, and tuning.
- [`vignette("algorithm", package = "somalign")`](https://mdmanurung.github.io/somalign/articles/algorithm.md):
  algorithm and output-column interpretation.

# somalign

`somalign` aligns query self-organising maps to fixed `kohonen`
reference maps using codebook-level unbalanced entropic optimal
transport.

The main label-transfer product is a transferred label, confidence,
margin, and acceptance flag for each query cell. Direct projection into
the reference SOM is the conservative reference-node assignment and the
preferred basis for composition/abundance summaries. Transport-corrected
projection columns are auxiliary diagnostics for visualisation,
annotation, and triage.

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
| `old_som_unit` | Direct reference-node assignment |
| `old_som_distance` | Distance to assigned reference node |
| `outside_reference_distance` | `TRUE` if distance exceeds reference quantile threshold |
| `outside_reference_surprisal` | Chi-squared-style anomaly score for the assigned reference node |
| `outside_reference_pvalue` | Approximate p-value for the anomaly score |
| `outside_reference_top_marker` | Marker contributing most to the anomaly score |
| `final_status` | `"inside_reference"` / `"outside_reference"` / `"unknown_reference_distance"` |
| `old_som_label` | Majority label of the assigned reference node |
| `old_som_label_confidence` | Fraction of that node’s training mass with the majority label |
| `corrected_som_unit` | OT-corrected reference-node assignment (auxiliary) |
| `correction_norm` | Magnitude of the OT correction shift |
| `transferred_label` | Top label from `correspondence %*% reference$label_prob` |
| `transferred_label_confidence` | Confidence of the transferred label |
| `transferred_label_second` | Runner-up transferred label |
| `transferred_label_margin` | Top minus runner-up transferred-label confidence |
| `transferred_label_accepted` | `TRUE` if match fraction and confidence thresholds are met |

## Learn more

- [`vignette("somalign", package = "somalign")`](https://mdmanurung.github.io/somalign/articles/somalign.md):
  quick start from reference training to result inspection.
- [`vignette("pretrained-old-and-new-soms", package = "somalign")`](https://mdmanurung.github.io/somalign/articles/pretrained-old-and-new-soms.md):
  existing SOMs, saved codebooks, diagnostics, and tuning.
- [`vignette("algorithm", package = "somalign")`](https://mdmanurung.github.io/somalign/articles/algorithm.md):
  algorithm and output-column interpretation.
- [`vignette("validating-label-transfer", package = "somalign")`](https://mdmanurung.github.io/somalign/articles/validating-label-transfer.md):
  label-transfer validation, calibration, tuning, and soft abundance
  summaries.

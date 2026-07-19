# Label-transfer accuracy metrics

Computes overall accuracy, macro-averaged F1, multiclass Matthews
correlation coefficient (Gorodkin's \\R_K\\), per-class
precision/recall/F1, the confusion matrix, and coverage (the fraction of
cells with an accepted prediction), given predicted and ground-truth
labels.

## Usage

``` r
somalign_label_metrics(predicted, truth, accepted = NULL)
```

## Arguments

- predicted:

  Character vector of predicted labels (`NA` = abstain).

- truth:

  Character vector of ground-truth labels, same length.

- accepted:

  Optional logical vector, same length. When supplied, only `accepted`
  (and non-`NA`) predictions are scored; the rest are abstentions. When
  `NULL` (default), all non-`NA` predictions are scored.

## Value

A list of class `somalign_label_metrics` with `accuracy`, `macro_f1`,
`mcc`, `per_class` (data frame), `confusion` (table), `n` (scored
predictions), `coverage`, and `accuracy_all`.

## Details

Metrics are computed on the *accepted* predictions only; `coverage`
reports what fraction of cells that is. `accuracy_all` additionally
scores abstentions (rejected or `NA` predictions) as wrong, for a
coverage-penalised view.

## Examples

``` r
truth <- rep(c("A", "B", "C"), each = 10)
pred  <- truth; pred[c(1, 12, 25)] <- "B"
somalign_label_metrics(pred, truth)
#> <somalign_label_metrics>
#>   accuracy = 0.9333  macro_f1 = 0.9346  MCC = 0.9045
#>   scored = 30  coverage = 100.0%  accuracy_all = 0.9333
```

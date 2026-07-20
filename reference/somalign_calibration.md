# Confidence calibration of label transfer

Bins predictions by a confidence-like score in \[0, 1\] (e.g. the
transfer confidence or margin) and compares mean score to empirical
accuracy per bin, yielding a reliability table plus the expected (ECE)
and maximum (MCE) calibration error and the top-label Brier score. A
well-calibrated model has mean score ~ accuracy in every bin (ECE ~ 0).

## Usage

``` r
somalign_calibration(score, correct, n_bins = 10L)
```

## Arguments

- score:

  Numeric vector in \[0, 1\]: the confidence/margin per prediction.

- correct:

  Logical vector, same length: whether each prediction was right.

- n_bins:

  Positive integer. Number of equal-width bins over \[0, 1\]. Default
  `10L`.

## Value

A list of class `somalign_calibration` with `table` (per-bin
`score_mean`, `accuracy`, `n`), `ece`, `mce`, `brier`, `n` (predictions
actually scored, after dropping `NA` score/correct pairs), `n_total`
(all supplied predictions), and `coverage` (`n / n_total`). Abstentions
(whose score/correct are `NA`) are dropped, so `ece`/`mce`/`brier`
describe the scored subset only; compare `coverage` before comparing
calibration across methods that abstain at different rates.

## Examples

``` r
set.seed(1)
score <- runif(200)
correct <- runif(200) < score        # perfectly calibrated by construction
somalign_calibration(score, correct)
#> <somalign_calibration>
#>   ECE = 0.0693  MCE = 0.1270  Brier = 0.1733  (scored n = 200, coverage = 100.0%)
#>   reliability (score_mean -> accuracy, n):
#>     0.06 -> 0.00  (12)
#>     0.15 -> 0.06  (17)
#>     0.25 -> 0.36  (22)
#>     0.36 -> 0.43  (23)
#>     0.46 -> 0.58  (24)
#>     0.55 -> 0.58  (19)
#>     0.65 -> 0.57  (21)
#>     0.75 -> 0.78  (27)
#>     0.86 -> 0.90  (20)
#>     0.95 -> 1.00  (15)
```

# Plot query vs reference SOM code ranges per marker

Visualises the min-to-max range of each marker's SOM codes for both the
query and reference codebooks. Colours indicate the overlap flag
computed by
[`somalign_check_codebook_alignment()`](https://mdmanurung.github.io/somalign/reference/somalign_check_codebook_alignment.md):
`ok` (green), `warning` (orange), or `critical` (red).

## Usage

``` r
somalign_plot_codebook_range(check)
```

## Arguments

- check:

  A `somalign_codebook_check` object, as returned by
  [`somalign_check_codebook_alignment()`](https://mdmanurung.github.io/somalign/reference/somalign_check_codebook_alignment.md).

## Value

A `ggplot` object.

## Examples

``` r
set.seed(1)
mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
              dimnames = list(NULL, c("F1", "F2")))
ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
                                rlen = 5)
qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
                      rlen = 5)
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
chk <- somalign_check_codebook_alignment(qry$codebook, ref,
                                         stop_if_critical = FALSE)
somalign_plot_codebook_range(chk)
```

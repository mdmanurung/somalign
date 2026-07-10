# Aligning SOMs with somalign

`somalign` keeps direct projection into the fixed reference SOM as the
primary classification. Optimal transport is used to describe
query-to-reference node correspondence and to produce auxiliary
corrected projection columns.

## Basic workflow

``` r

library(kohonen)
library(somalign)

set.seed(1)
old <- rbind(
  matrix(rnorm(80, mean = -1), ncol = 4),
  matrix(rnorm(80, mean = 1), ncol = 4)
)
colnames(old) <- paste0("f", seq_len(ncol(old)))
labels <- rep(c("low", "high"), each = 20)

reference <- somalign_train_reference(
  old,
  labels = labels,
  grid = kohonen::somgrid(2, 2, "hexagonal"),
  rlen = 5
)

query <- old + 0.2
query_obj <- somalign_query(
  query,
  reference,
  grid = kohonen::somgrid(2, 2, "hexagonal"),
  rlen = 5
)

fit <- somalign_fit(query_obj, reference, solver = "internal")
results <- somalign_results(fit)
head(results)
#>   sample_id query_som_unit old_som_unit old_som_distance
#> 1         1              4            1        1.6019468
#> 2         2              1            4        0.8377855
#> 3         3              1            4        1.3090895
#> 4         4              3            4        1.9063036
#> 5         5              1            4        0.8813879
#> 6         6              1            4        0.9746628
#>   old_som_distance_threshold outside_reference_distance     final_status
#> 1                   1.775302                      FALSE inside_reference
#> 2                   1.919417                      FALSE inside_reference
#> 3                   1.919417                      FALSE inside_reference
#> 4                   1.919417                      FALSE inside_reference
#> 5                   1.919417                      FALSE inside_reference
#> 6                   1.919417                      FALSE inside_reference
#>   old_som_label old_som_label_confidence corrected_som_unit
#> 1          high                      0.6                  1
#> 2           low                      1.0                  4
#> 3           low                      1.0                  4
#> 4           low                      1.0                  4
#> 5           low                      1.0                  4
#> 6           low                      1.0                  4
#>   corrected_som_distance corrected_som_distance_threshold
#> 1              0.2975369                         1.775302
#> 2              0.6852730                         1.919417
#> 3              1.2155549                         1.919417
#> 4              1.7911076                         1.919417
#> 5              0.8350028                         1.919417
#> 6              0.5477944                         1.919417
#>   corrected_outside_reference_distance correction_norm transferred_label
#> 1                                FALSE       1.3557100              high
#> 2                                FALSE       0.4712923               low
#> 3                                FALSE       0.4712923               low
#> 4                                FALSE       1.1474817              high
#> 5                                FALSE       0.4712923               low
#> 6                                FALSE       0.4712923               low
#>   transferred_label_confidence transferred_label_accepted
#> 1                    0.6000000                       TRUE
#> 2                    1.0000000                       TRUE
#> 3                    1.0000000                       TRUE
#> 4                    0.9166529                       TRUE
#> 5                    1.0000000                       TRUE
#> 6                    1.0000000                       TRUE
```

## Output columns

[`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md)
returns one row per query sample. The columns fall into three groups:

**Direct projection** (primary result; use these for downstream
analysis):

- `old_som_unit`: reference node assigned by direct projection
- `old_som_distance`: Euclidean distance to that node in
  reference-scaled space
- `outside_reference_distance`: `TRUE` if distance exceeds the per-node
  distance quantile threshold derived from reference training data
- `final_status`: `"inside_reference"`, `"outside_reference"`, or
  `"unknown_reference_distance"` (when no threshold is available for the
  node)
- `old_som_label`: majority label of the assigned reference node
- `old_som_label_confidence`: fraction of that node’s training mass with
  the majority label

**OT-corrected projection** (auxiliary; use for visualisation and
triage):

- `corrected_som_unit`: reference node after applying the OT correction
  shift
- `corrected_som_distance`: distance to the corrected node
- `correction_norm`: magnitude of the correction vector in
  reference-scaled space

**Label transfer** (auxiliary, based on OT correspondence):

- `transferred_label`: label transferred via OT node correspondence
- `transferred_label_confidence`: confidence of that transferred label
- `transferred_label_accepted`: `TRUE` when match fraction ≥
  `min_match_fraction` and label confidence ≥ `confidence_threshold`

``` r

names(results)
#>  [1] "sample_id"                           
#>  [2] "query_som_unit"                      
#>  [3] "old_som_unit"                        
#>  [4] "old_som_distance"                    
#>  [5] "old_som_distance_threshold"          
#>  [6] "outside_reference_distance"          
#>  [7] "final_status"                        
#>  [8] "old_som_label"                       
#>  [9] "old_som_label_confidence"            
#> [10] "corrected_som_unit"                  
#> [11] "corrected_som_distance"              
#> [12] "corrected_som_distance_threshold"    
#> [13] "corrected_outside_reference_distance"
#> [14] "correction_norm"                     
#> [15] "transferred_label"                   
#> [16] "transferred_label_confidence"        
#> [17] "transferred_label_accepted"
```

## Diagnostics

``` r

diag <- somalign_diagnostics(fit)

# Which solver was used, and any fallback notes
diag$solver
#> $requested
#> [1] "internal"
#> 
#> $used
#> [1] "internal"
#> 
#> $notes
#> character(0)
#> 
#> $iterations
#> [1] 164
#> 
#> $epsilon
#> [1] 0.05
#> 
#> $rho_query
#> [1] 1
#> 
#> $rho_ref
#> [1] 1

# Per-node OT statistics
diag$ot
#> $transport_mass
#> [1] 0.7036713
#> 
#> $row_mass
#>         V1         V2         V3         V4 
#> 0.29095033 0.27598507 0.09991421 0.03682170 
#> 
#> $col_mass
#>         V1         V2         V3         V4 
#> 0.05160523 0.17396174 0.18715271 0.29095163 
#> 
#> $query_mass
#> [1] 0.300 0.475 0.175 0.050
#> 
#> $reference_mass
#> [1] 0.125 0.150 0.300 0.425
#> 
#> $match_fraction
#> [1] 0.9698344 0.5810212 0.5709383 0.7364339
#> 
#> $match_mass_ratio
#> [1] 0.9698344 0.5810212 0.5709383 0.7364339
#> 
#> $max_row_mass_error
#> [1] 0.1990149
#> 
#> $max_col_mass_error
#> [1] 0.1340484
```

## Sensitivity analysis

Use
[`somalign_sensitivity_grid()`](https://mdmanurung.github.io/somalign/reference/somalign_sensitivity_grid.md)
to evaluate robustness across OT hyperparameter combinations before
committing to a single fit:

``` r

grid_results <- somalign_sensitivity_grid(
  query_obj,
  reference,
  epsilon   = c(0.05, 0.1, 0.5),
  rho_query = c(0.5, 1),
  rho_ref   = 1,
  solver    = "internal"
)
#> somalign_fit: 3 query node(s) have match_mass_ratio > 1 (max 1.42); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 3 query node(s) have match_mass_ratio > 1 (max 1.27); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
grid_results
#>   epsilon rho_query rho_ref   solver transport_mass mean_match_fraction
#> 1    0.05       0.5       1 internal      0.6366066           0.6611138
#> 2    0.10       0.5       1 internal      0.6755324           0.6963160
#> 3    0.50       0.5       1 internal      1.0182705           0.9632075
#> 4    0.05       1.0       1 internal      0.7036713           0.7145570
#> 5    0.10       1.0       1 internal      0.7352817           0.7450841
#> 6    0.50       1.0       1 internal      1.0123260           0.9727605
#>   max_row_mass_error max_col_mass_error accepted_label_fraction
#> 1         0.24518349         0.14251734                    1.00
#> 2         0.22533288         0.12950682                    1.00
#> 3         0.06990570         0.06786276                    0.75
#> 4         0.19901493         0.13404837                    1.00
#> 5         0.18270660         0.12510911                    1.00
#> 6         0.05175497         0.07290494                    0.75
#>   outside_direct_fraction outside_corrected_fraction
#> 1                   0.175                      0.100
#> 2                   0.175                      0.100
#> 3                   0.175                      0.050
#> 4                   0.175                      0.125
#> 5                   0.175                      0.125
#> 6                   0.175                      0.050
```

Each row reports `transport_mass`, `mean_match_fraction`, and marginal
mass errors for one parameter combination. Consistent values across rows
indicate the result is not sensitive to the choice of `epsilon` or
`rho`.

## OT hyperparameters

[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
exposes three OT hyperparameters:

- **`epsilon`** (default `0.05`): entropic regularisation. Smaller
  values give a sparser, more deterministic transport plan; larger
  values give a smoother plan that converges faster. If \> 1% of
  Sinkhorn kernel entries underflow, a warning gives a safe lower bound.
- **`rho_query`** / **`rho_ref`** (default `1`): KL-divergence marginal
  penalty. Larger values enforce tighter mass conservation (approaches
  balanced OT); smaller values allow more mass to be discarded
  (appropriate when the query contains populations absent from the
  reference). The deviation of row sums from input masses is reported in
  `diagnostics$ot$max_row_mass_error`.

## Large datasets

For large query sets, projection is computed in chunks to cap peak
memory:

``` r

fit <- somalign_fit(query_obj, reference, chunk_size = 5000L)
```

The default `chunk_size = 10000L` limits peak allocation to
`10000 × n_reference_nodes` doubles. Set `chunk_size = Inf` to use a
single full matrix allocation.

# Package index

## Building reference and query objects

- [`somalign_train_reference()`](https://mdmanurung.github.io/somalign/reference/somalign_train_reference.md)
  : Train a reference SOM and build a somalign reference
- [`somalign_reference()`](https://mdmanurung.github.io/somalign/reference/somalign_reference.md)
  : Build a reference object from an existing SOM and old data
- [`somalign_reference_from_nodes()`](https://mdmanurung.github.io/somalign/reference/somalign_reference_from_nodes.md)
  : Build a reference object from saved node-level artifacts
- [`somalign_reference_from_som()`](https://mdmanurung.github.io/somalign/reference/somalign_reference_from_som.md)
  : Build a reference object directly from a trained kohonen SOM
- [`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md)
  : Prepare query data and attach or train a query SOM
- [`somalign_query_from_som()`](https://mdmanurung.github.io/somalign/reference/somalign_query_from_som.md)
  : Build a query object from a pre-trained kohonen SOM

## Preprocessing

- [`somalign_normalize()`](https://mdmanurung.github.io/somalign/reference/somalign_normalize.md)
  : Global pre-correction of query data to match the reference
  distribution
- [`somalign_quantile_normalize()`](https://mdmanurung.github.io/somalign/reference/somalign_quantile_normalize.md)
  : Divide each feature column by its upper quantile

## Fitting

- [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
  : Align a query SOM to a reference SOM
- [`somalign_fit_anchored()`](https://mdmanurung.github.io/somalign/reference/somalign_fit_anchored.md)
  : Align a query SOM to a reference SOM using anchor sample pairs
- [`somalign_fit_two_pass()`](https://mdmanurung.github.io/somalign/reference/somalign_fit_two_pass.md)
  : Two-pass alignment decomposing correction into global and local
  components

## Results

- [`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md)
  : Return per-sample somalign results

## Diagnostics and sensitivity

- [`somalign_diagnostics()`](https://mdmanurung.github.io/somalign/reference/somalign_diagnostics.md)
  : Extract somalign diagnostics
- [`somalign_sensitivity_grid()`](https://mdmanurung.github.io/somalign/reference/somalign_sensitivity_grid.md)
  : Run an OT sensitivity grid
- [`somalign_som_stability()`](https://mdmanurung.github.io/somalign/reference/somalign_som_stability.md)
  : Assess alignment stability across query SOM random seeds
- [`somalign_worst_nodes()`](https://mdmanurung.github.io/somalign/reference/somalign_worst_nodes.md)
  : Return worst-projecting query SOM nodes
- [`somalign_check_codebook_alignment()`](https://mdmanurung.github.io/somalign/reference/somalign_check_codebook_alignment.md)
  : Check query-reference codebook alignment before cell-level
  computation

## Plots

- [`somalign_plot_codebook_range()`](https://mdmanurung.github.io/somalign/reference/somalign_plot_codebook_range.md)
  : Plot query vs reference SOM code ranges per marker
- [`somalign_plot_correction()`](https://mdmanurung.github.io/somalign/reference/somalign_plot_correction.md)
  : Plot per-node correction norms
- [`somalign_plot_label_confusion()`](https://mdmanurung.github.io/somalign/reference/somalign_plot_label_confusion.md)
  : Plot label transfer confusion heatmap
- [`somalign_plot_marker_distributions()`](https://mdmanurung.github.io/somalign/reference/somalign_plot_marker_distributions.md)
  : Plot per-marker cell distributions before projection
- [`somalign_plot_mass_balance()`](https://mdmanurung.github.io/somalign/reference/somalign_plot_mass_balance.md)
  : Plot node mass balance
- [`somalign_plot_match_fraction()`](https://mdmanurung.github.io/somalign/reference/somalign_plot_match_fraction.md)
  : Plot per-node match fraction
- [`somalign_plot_outside_fraction()`](https://mdmanurung.github.io/somalign/reference/somalign_plot_outside_fraction.md)
  : Plot fraction of cells outside reference thresholds

## Print methods

- [`print(`*`<somalign_anchored_fit>`*`)`](https://mdmanurung.github.io/somalign/reference/print.somalign_anchored_fit.md)
  : Print a somalign_anchored_fit object
- [`print(`*`<somalign_fit>`*`)`](https://mdmanurung.github.io/somalign/reference/print.somalign_fit.md)
  : Print a somalign_fit object
- [`print(`*`<somalign_query>`*`)`](https://mdmanurung.github.io/somalign/reference/print.somalign_query.md)
  : Print a somalign_query object
- [`print(`*`<somalign_reference>`*`)`](https://mdmanurung.github.io/somalign/reference/print.somalign_reference.md)
  : Print a somalign_reference object

# Package index

## Training and building inputs

- [`somalign_train_reference()`](https://mdmanurung.github.io/somalign/reference/somalign_train_reference.md)
  : Train a reference SOM and build a somalign reference
- [`somalign_reference()`](https://mdmanurung.github.io/somalign/reference/somalign_reference.md)
  : Build a reference object from an existing SOM and old data
- [`somalign_reference_from_nodes()`](https://mdmanurung.github.io/somalign/reference/somalign_reference_from_nodes.md)
  : Build a reference object from saved node-level artifacts
- [`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md)
  : Prepare query data and attach or train a query SOM

## Fitting and results

- [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
  : Align a query SOM to a reference SOM
- [`somalign_fit_anchored()`](https://mdmanurung.github.io/somalign/reference/somalign_fit_anchored.md)
  : Align a query SOM to a reference SOM using anchor sample pairs
- [`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md)
  : Return per-sample somalign results

## Diagnostics and sensitivity

- [`somalign_diagnostics()`](https://mdmanurung.github.io/somalign/reference/somalign_diagnostics.md)
  : Extract somalign diagnostics
- [`somalign_sensitivity_grid()`](https://mdmanurung.github.io/somalign/reference/somalign_sensitivity_grid.md)
  : Run an OT sensitivity grid
- [`somalign_som_stability()`](https://mdmanurung.github.io/somalign/reference/somalign_som_stability.md)
  : Assess alignment stability across query SOM random seeds

## Print methods

- [`print(`*`<somalign_anchored_fit>`*`)`](https://mdmanurung.github.io/somalign/reference/print.somalign_anchored_fit.md)
  : Print a somalign_anchored_fit object
- [`print(`*`<somalign_fit>`*`)`](https://mdmanurung.github.io/somalign/reference/print.somalign_fit.md)
  : Print a somalign_fit object
- [`print(`*`<somalign_query>`*`)`](https://mdmanurung.github.io/somalign/reference/print.somalign_query.md)
  : Print a somalign_query object
- [`print(`*`<somalign_reference>`*`)`](https://mdmanurung.github.io/somalign/reference/print.somalign_reference.md)
  : Print a somalign_reference object

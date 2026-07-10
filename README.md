# somalign

`somalign` aligns query self-organising maps to fixed `kohonen` reference maps
with codebook-level unbalanced entropic optimal transport.

The conservative primary result is direct projection into the old reference SOM.
Transport-corrected projections are returned as auxiliary columns for
visualisation, annotation, and triage.

```r
library(kohonen)
library(somalign)

reference <- somalign_train_reference(old_matrix, labels = old_labels)
query <- somalign_query(new_matrix, reference)
fit <- somalign_fit(query, reference)
results <- somalign_results(fit)
```

Python POT is optional. When `reticulate` can import `ot.unbalanced`,
`solver = "auto"` uses POT; otherwise `somalign` falls back to an internal
generalized Sinkhorn solver and records that choice in diagnostics.

If you already trained a new/query SOM, train it on new samples transformed
with the old reference `center` and `scale`, then pass it as
`somalign_query(new_data, reference, som_query = new_som)`. The canonical old
projection is returned in `old_som_unit`; the auxiliary corrected old-node
projection is returned in `corrected_som_unit`.

If you build a reference from an existing old SOM, state the codebook
coordinate system explicitly:

```r
reference <- somalign_reference(
  old_som,
  old_matrix,
  codebook_space = "reference_scaled" # or "raw"
)
```

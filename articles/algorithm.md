# Algorithm and interpretation

`somalign` aligns a query SOM to a fixed reference SOM using
codebook-level KL-unbalanced entropic optimal transport. This page walks
through each stage of the pipeline, explaining how direct projection,
the OT correspondence, correction vectors, label transfer, and soft
abundance estimates each contribute to the current API.

The OT stage is implemented inside the package. `solver = "internal"` is
the default primal-domain Sinkhorn iteration, while
`solver = "log_domain"` and `solver = "annealing"` provide more stable
log-potential paths for small `epsilon`, label-guided fits, or difficult
cost matrices.

The primary outputs depend on the task. For query-cell annotation, use
the transferred label, confidence, margin, and acceptance flag from
stage 5. For reference-node assignment and cross-batch composition, use
direct projection from stage 1. For abundance profiles, prefer the soft
frequency estimator from stage 6. The corrected projection columns from
stages 3 and 4 are auxiliary diagnostics, shown muted in the diagram
below.

> **Correction vectors are a diagnostic path.** The barycentric
> correction (stages 3 and 4 below) can over-merge distinct populations
> at larger epsilon: a conditional-expectation contraction that
> [`somalign_topology_audit()`](https://mdmanurung.github.io/somalign/reference/somalign_topology_audit.md)
> quantifies. Use the corrected projection columns for triage and
> visualisation. If you need a cell-by-marker corrected expression
> matrix, use the dedicated
> [`somalign_correct_expression()`](https://mdmanurung.github.io/somalign/reference/somalign_correct_expression.md)
> path on an anchored subspace-restricted fit or a two-pass fit, and
> audit topology first. Label transfer (stage 5) is computed from the
> transport correspondence and never uses corrected coordinates.

## Pipeline overview

------------------------------------------------------------------------

## Stage 0: Training both SOMs

The reference SOM is trained on labelled old data after z-score
normalisation with the old-data mean $`\mu`$ and standard deviation
$`\sigma`$. The query SOM must be trained on new data transformed with
those same reference parameters (not new-data z-scores), so that both
codebooks lie in the same feature coordinate system.

![](algorithm_files/figure-html/stage0-1.png)

------------------------------------------------------------------------

## Stage 1: Direct projection (reference-node assignment)

Every new sample $`x_s`$ is assigned to its nearest reference node by
Euclidean distance:
``` math
k^* = \arg\min_k \|x_s - C_{\text{ref},k}\|
```

No transport is involved. The assignment is deterministic given the
reference codebook. Use it when the question is “which reference node is
this cell nearest to?” and for composition or abundance summaries
(`old_som_unit` / `old_som_label`). Label transfer is handled separately
in stage 5.

![](algorithm_files/figure-html/stage1-1.png)

**Output:** `old_som_unit = k*`, `old_som_distance`,
`old_som_distance_threshold`, `outside_reference_distance`,
`outside_reference_surprisal`, `outside_reference_pvalue`,
`outside_reference_top_marker`, `final_status`, `old_som_label`,
`old_som_label_confidence`.

------------------------------------------------------------------------

## Stage 2: OT on codebooks

Transport is solved once on the **codebooks** ($`M`$ query nodes × $`K`$
reference nodes), not on individual samples, keeping the problem
dimensionality independent of sample count.

### 2a: Cost matrix

![](algorithm_files/figure-html/stage2a-1.png)

### 2b: OT objective

The transport plan $`P`$ minimises:

``` math
\sum_{ij} P_{ij} C_{ij}
  + \varepsilon \cdot \mathrm{KL}(P \,\|\, \mathbf{1})
  + \rho_q \cdot \mathrm{KL}(P\mathbf{1} \,\|\, \mathbf{a})
  + \rho_r \cdot \mathrm{KL}(P^\top\!\mathbf{1} \,\|\, \mathbf{b})
```

where $`\mathbf{a}`$ and $`\mathbf{b}`$ are the query and reference node
masses. Before the solve, the squared-distance cost matrix is divided by
its median positive entry, so `epsilon` is on a roughly stable scale
across marker panels. Optional feature weights rescale marker columns
before computing the cost. Optional diagonal boost, anchor cost bonuses,
and label-guided penalties then modify the normalised cost.

The entropic term $`\varepsilon \cdot \mathrm{KL}(P \,\|\, \mathbf{1})`$
yields the Gibbs kernel $`K_{ij} = \exp(-C_{ij}/\varepsilon)`$ used by
the primal-domain solver. The log-domain and annealing solvers use
equivalent log-potential updates and avoid explicitly materialising an
underflowing kernel. The KL marginal penalties relax the balance
constraint: mass can be destroyed rather than transported, so a query
node with no good reference counterpart routes only a small fraction of
its mass rather than forcing an artificial match. The plan is found by
the internal pure-R Sinkhorn scaling implementation selected by
`solver`.

### 2c: Transport plan

![](algorithm_files/figure-html/stage2c-1.png)

`match_fraction` near 1 means the node transported most of its mass to
reference nodes; near 0 means little overlap with any reference
population. Label transfer is rejected for low-match nodes.

------------------------------------------------------------------------

## Stage 3: Correction vectors

For each query node $`i`$ the OT plan defines a weighted mean
displacement toward its matched reference nodes:

``` math
\Delta_i = \frac{\displaystyle\sum_j P_{ij}\,(C_{\text{ref},j} - C_{\text{query},i})}
                  {\displaystyle\sum_j P_{ij}}
```

![](algorithm_files/figure-html/stage3-1.png)

A small $`\|\Delta_i\|`$ means the query node sits close to its matched
reference positions. A large value indicates a systematic offset
(potentially a batch effect or a true population shift) and is worth
examining before accepting the corrected projection.

------------------------------------------------------------------------

## Stage 4: Corrected projection (auxiliary)

Every sample in query node $`i`$ is shifted by the same $`\Delta_i`$ and
re-projected:

``` math
\text{corrected\_som\_unit}_s = \arg\min_k \|(x_s + \Delta_i) - C_{\text{ref},k}\|
```

![](algorithm_files/figure-html/stage4-1.png)

Both assignments are reported in the results table. `old_som_unit`
(blue) is the conservative direct reference-node assignment;
`corrected_som_unit` (red), `corrected_som_distance`, and
`corrected_outside_reference_distance` are intended for visualisation
and triage, not as replacement assignments.

------------------------------------------------------------------------

## Stage 5: Label transfer (query-cell annotation)

For each query node $`i`$, the raw transport row is first normalised
over reference nodes:

``` math
T_{ij} = \frac{P_{ij}}{\sum_l P_{il}}.
```

The query-node label posterior is then the full transported mixture of
reference-node label probabilities:

``` math
\pi_i = T_i L_{\text{ref}},
```

where $`L_{\text{ref}}[j,c] = \Pr(c \mid \text{reference node } j)`$.
The transferred label is `argmax_c pi_i[c]`; confidence is the top
posterior probability; the second label and margin expose close calls.
This is more accurate than treating the single highest-mass reference
node as the only donor.

![](algorithm_files/figure-html/stage5-1.png)

Transfer is rejected at two gates: `match_fraction < 0.05` (node has no
meaningful reference overlap) and `transferred_label_confidence < 0.60`
(the transported label posterior is too mixed to assign a clean label).
The result table also exposes `transferred_label_second`,
`transferred_label_second_confidence`, and `transferred_label_margin` so
borderline cells can be triaged without discarding the posterior
information.

------------------------------------------------------------------------

## Stage 6: Soft label and frequency projection

[`somalign_soft_labels()`](https://mdmanurung.github.io/somalign/reference/somalign_soft_labels.md)
and
[`somalign_soft_frequencies()`](https://mdmanurung.github.io/somalign/reference/somalign_soft_frequencies.md)
are separate from OT label transfer. They project each query cell
directly to the reference codebook with a Gaussian kernel over its `k`
nearest reference nodes, then average the corresponding reference label
probabilities or user-supplied node-group indicators. No acceptance gate
or out-of-reference gate is applied: every query cell contributes to the
soft distribution unless all of its neighbouring reference nodes are
unlabelled.

This is the soft analogue of hard direct projection. It is mainly useful
for per-sample cluster or metacluster abundance, where hard nearest-node
counts are sensitive to small boundary shifts. Soft projection changes
the frequency estimate; it does not change the most-likely hard label in
[`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md).

``` r

soft_labels <- somalign_soft_labels(fit, k = 8)

# sample_id: one group label per query cell
soft_freq <- somalign_soft_frequencies(fit, group = sample_id, k = 8)

# node2metacluster: one metacluster label per reference node
soft_meta <- somalign_soft_frequencies(
  fit,
  group = sample_id,
  node_groups = node2metacluster,
  normalize = FALSE
)
```

------------------------------------------------------------------------

## Limitations of the OT correction and uncertainty scores

Three structural limitations govern when corrected projection columns
and uncertainty scores should and should not be trusted. These
limitations do not change the direct projection, and label transfer
should be interpreted through its own confidence, margin, and
calibration diagnostics.

### Barycentric contraction (F1)

The correction target for node $`i`$ is the transport-weighted mean of
matched reference codebook vectors:

``` math
\hat{C}_i = \frac{\sum_j P_{ij}\,C_{\text{ref},j}}{\sum_j P_{ij}}
```

This barycenter is a **conditional mean**, so the corrected positions
always lie inside the convex hull of the reference codebook. Two
consequences follow. First, variance shrinks: query nodes that genuinely
sit near the periphery of a reference cluster will be pulled inward.
Second, dense reference regions attract mass. For example, if the
reference has one very large T-cell cluster and one small NK cluster,
even a predominantly NK query node will have its correction vector
biased toward the T-cell centroid. Contraction grows as `epsilon`
increases, because a larger entropic penalty spreads each transport row
$`P[i,
\cdot]`$ across more reference nodes. Interpret `corrected_som_unit`
with `correction_norm`: a large shift ($`\|\Delta_i\|`$) in a
high-epsilon fit is the most suspect combination.

### Unpaired OT without anchors (F2)

OT operates on **node mass distributions**, not on individual samples.
If the query batch consists of literally remeasured biological units
(the same cells or patients measured again), that ground-truth
sample-level correspondence is available but unused here. Anchor-based
methods such as CytoNorm exploit exactly this pairing, using anchor
controls measured in both old and new batches to identify the batch map
with less ambiguity. The plain
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
path is designed for the common case where no such anchors exist: only
the phenotypic structure of the populations, encoded in the codebooks,
is shared. If you do have remeasured anchor controls, use
[`somalign_fit_anchored()`](https://mdmanurung.github.io/somalign/reference/somalign_fit_anchored.md).
The `correction = "cost_bonus"` mode biases the OT cost toward
anchor-supported routes; `correction = "subspace"` and `"both"` estimate
an anchor-derived batch subspace and restrict correction vectors to that
subspace.

### Novelty and confidence are parameter-dependent (F3)

`match_fraction` (and the derived label-transfer gate) reflects the
$`\varepsilon`$/$`\rho`$ regime as much as any biological novelty. A low
`match_fraction` means little query mass reached reference nodes, but
that can happen because (a) the query contains a genuinely novel
population, or (b) `rho_query` is small, (c) `epsilon` is large, or (d)
the query SOM topology differs from the reference. There is no null
model calibrating how much destroyed mass constitutes evidence of
novelty.

`outside_reference_surprisal` and `outside_reference_pvalue` provide a
separate per-cell anomaly score under a diagonal-Gaussian reference-node
model. That p-value is useful for ranking unusual cells and flagging
large marker-level deviations, but it is approximate and can be
anti-conservative for heavy-tailed marker distributions. Use
[`somalign_sensitivity_grid()`](https://mdmanurung.github.io/somalign/reference/somalign_sensitivity_grid.md),
[`somalign_epsilon_sweep()`](https://mdmanurung.github.io/somalign/reference/somalign_epsilon_sweep.md),
and
[`somalign_calibration()`](https://mdmanurung.github.io/somalign/reference/somalign_calibration.md)
to assess whether low match, high surprisal, and low confidence are
stable before interpreting them as novel populations.

------------------------------------------------------------------------

## Output columns

The results table contains direct projection columns, transferred-label
columns, and auxiliary corrected-projection columns. Use the direct
columns for reference-node assignment and abundance, the
transferred-label columns for query-cell annotation, and the corrected
columns only as correction diagnostics alongside
[`somalign_diagnostics()`](https://mdmanurung.github.io/somalign/reference/somalign_diagnostics.md).

| Column | Stage | Purpose |
|----|----|----|
| `sample_id` | 0 | Query sample identifier |
| `query_som_unit` | 0 | Query SOM node containing the sample |
| `old_som_unit` | 1 | Primary node assignment |
| `old_som_distance` | 1 | Distance to that node |
| `old_som_distance_threshold` | 1 | Direct distance threshold |
| `outside_reference_distance` | 1 | Novelty flag |
| `outside_reference_surprisal` | 1 | Diagonal-Gaussian anomaly score for the assigned reference node |
| `outside_reference_pvalue` | 1 | Approximate p-value for the anomaly score |
| `outside_reference_pvalue_flag` | 1 | Optional p-value flag when `outside_pvalue_threshold` is supplied |
| `outside_reference_top_marker` | 1 | Marker contributing most to the anomaly score |
| `final_status` | 1 | inside / outside / unknown |
| `old_som_label` | 1 | Primary label |
| `old_som_label_confidence` | 1 | Label purity at node |
| `corrected_som_unit` | 4 | OT-corrected assignment |
| `corrected_som_distance` | 4 | Distance after applying the node correction |
| `corrected_som_distance_threshold` | 4 | Corrected distance threshold |
| `corrected_outside_reference_distance` | 4 | Novelty flag after correction |
| `correction_norm` | 3/4 | Shift magnitude $`\|\Delta_i\|`$ |
| `transferred_label` | 5 | Top label from `correspondence %*% reference$label_prob` |
| `transferred_label_confidence` | 5 | Top transported-label posterior probability |
| `transferred_label_accepted` | 5 | Gate result (TRUE / FALSE) |
| `transferred_label_second` | 5 | Runner-up transported label |
| `transferred_label_second_confidence` | 5 | Runner-up transported-label posterior probability |
| `transferred_label_margin` | 5 | Top confidence minus runner-up confidence |

OT runs once on the $`M \times K`$ codebook matrix. The per-sample cost
is an $`O(n \cdot K)`$ nearest-node search, performed in configurable
chunks by
[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
to keep memory use predictable.

## Session info

    #> R version 4.6.1 (2026-06-24)
    #> Platform: x86_64-pc-linux-gnu
    #> Running under: Ubuntu 24.04.4 LTS
    #> 
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    #> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    #> 
    #> locale:
    #>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    #>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    #>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    #> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    #> 
    #> time zone: UTC
    #> tzcode source: system (glibc)
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices utils     datasets  methods   base     
    #> 
    #> other attached packages:
    #> [1] BiocStyle_2.40.0
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] cli_3.6.6           knitr_1.51          rlang_1.3.0        
    #>  [4] xfun_0.60           otel_0.2.0          jsonlite_2.0.0     
    #>  [7] glue_1.8.1          htmltools_0.5.9     sass_0.4.10        
    #> [10] rmarkdown_2.31      evaluate_1.0.5      jquerylib_0.1.4    
    #> [13] visNetwork_2.1.4    fastmap_1.2.0       yaml_2.3.12        
    #> [16] lifecycle_1.0.5     bookdown_0.47       BiocManager_1.30.27
    #> [19] DiagrammeR_1.0.12   compiler_4.6.1      fs_2.1.0           
    #> [22] RColorBrewer_1.1-3  htmlwidgets_1.6.4   digest_0.6.39      
    #> [25] R6_2.6.1            magrittr_2.0.5      bslib_0.11.0       
    #> [28] tools_4.6.1         pkgdown_2.2.1       cachem_1.1.0       
    #> [31] desc_1.4.3

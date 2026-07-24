# Plot achieved vs requested correlation heatmap

Side-by-side heatmaps comparing the requested and achieved correlation
matrices from a simulation object.

## Usage

``` r
plot_correlation_heatmap(sim, difference = FALSE)
```

## Arguments

- sim:

  An `aucmat_simulation` object.

- difference:

  If `TRUE`, show a third panel with the difference matrix (achieved -
  requested). Default `FALSE`.

## Value

A `ggplot2` object.

## Examples

``` r
# \donttest{
set.seed(42)
sim <- simulate_auc_matrix(n=50, prevalence=0.3,
  target_aucs=c(0.8,0.7,0.6), correlation=0.3, structure="exchangeable")
plot_correlation_heatmap(sim)

# }
```

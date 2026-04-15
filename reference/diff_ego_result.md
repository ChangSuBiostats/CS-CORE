# GO enrichment results from the differential co-expression vignette

GO enrichment results (from
[`clusterProfiler::enrichGO`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html))
for all WGCNA modules identified in the differential co-expression
analysis of CD14 Monocytes (COVID-19 vs. healthy controls). Saved for
rendering the vignette without re-running the full analysis.

## Usage

``` r
diff_ego_result
```

## Format

### `diff_ego_result`

A list of `enrichResult` objects, one per module.

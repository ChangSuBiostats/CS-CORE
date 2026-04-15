## Script to generate top_enrich_go_diff.rda
##
## Run this script once after completing the differential co-expression analysis
## in vignettes/differential_coexpression.Rmd to cache the GO enrichment results
## so that the vignette can render without re-running the full pipeline.
##
## Prerequisites: the objects `ego_result` and `top_enrich_clusters_diff` must
## be available in your R session (produced by the GO enrichment step in the vignette).

top_enrich_go_diff <- lapply(top_enrich_clusters_diff, function(i) ego_result[[i]]@result[1:3, ])

usethis::use_data(top_enrich_go_diff, overwrite = TRUE)

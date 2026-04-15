## Script to generate cached data objects for the differential co-expression vignette.
##
## Prerequisites: run the full analysis in vignettes/differential_coexpression.Rmd.orig
## first so that the following objects are available in your R session:
##   diff_result  -- output of diff_CSCORE()
##   ego_result   -- list of enrichResult objects from enrichGO()
##   module_list  -- list of gene vectors from WGCNA clustering

# module_list: gene module assignments from WGCNA
diff_module_list <- module_list

# ego_result: GO enrichment results for all modules
diff_ego_result <- ego_result

# diff_est_example: co-expression matrices for the example module (module 4),
# subsetted to the genes in that module to keep the object compact
gene_set <- module_list[[4]]
diff_est_example <- list(
  Healthy = diff_result$est_group1[gene_set, gene_set],
  COVID   = diff_result$est_group2[gene_set, gene_set]
)

usethis::use_data(diff_module_list, overwrite = TRUE)
usethis::use_data(diff_ego_result,  overwrite = TRUE)
usethis::use_data(diff_est_example, overwrite = TRUE)

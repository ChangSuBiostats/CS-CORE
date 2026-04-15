.get_assay_data <- function(object, assay = "RNA", layer = "counts") {
  # Try SeuratObject v5-style first
  out <- tryCatch(
    SeuratObject::GetAssayData(object = object, assay = assay, layer = layer),
    error = function(e) e
  )
  if (!inherits(out, "error")) return(out)

  # Fall back to v4-style slot
  SeuratObject::GetAssayData(object = object, assay = assay, slot = layer)
}

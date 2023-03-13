#' Function to be called to normalize single cell data
#'
#' @param seuratobject the seurat object you want normalized
#' @param method whichever normalization method you want to implement like "LogNormalize"
#' @param scale scale
#'
#' @return a seurat object with normalized count information
#' @export

normalize_sc <- function(seuratobject, method, scale) {
  Seurat::NormalizeData(seuratobject, normalization.method = method, scale.factor = scale)
}

#' Normalize an entire seurat list
#'
#' @param seuratlist your list of seurat objects you want normalized
#' @param method the normalization method of your choice
#' @param scale scale
#'
#' @return a list of seurat objects with normalized counts
#' @export

normalize_seurat_list <- function(seuratlist, method, scale){
  seulist.normalized <- lapply(seuratlist, normalize_sc, method=method, scale=scale)
  return(seulist.normalized)
}

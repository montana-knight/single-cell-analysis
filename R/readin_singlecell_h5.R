#' Import a single cell h5 file
#'
#' @param filename Your h5 file and any needed file path
#' @param sampleid What you want the file to be called upon import
#'
#' @return A seurat object
#' @export


readin_singlecell_h5 <- function(filename, sampleid) {
  Seurat::Read10X_h5(filename) %>%
    Seurat::CreateSeuratObject(
      project = sampleid) -> seu.obj
  return(seu.obj)
}

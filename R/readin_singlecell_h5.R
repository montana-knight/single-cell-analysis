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

#' Import multiple h5 files
#'
#' @param filenames The files you want imported, including necessary path information
#' @param samplenames What you want the samples to be called
#'
#' @return A list of seurat objects that correspond to your h5 files
#' @export
#'

readin_multiple_h5files <- function(filenames, samplenames){
  seulist <- list()
  for (i in 1:length(filenames)) {
    filename <- filenames[i]
    sampleid <- samplenames[i]
    seulist[[i]] <- readin_singlecell_h5(filename, sampleid)
    names(seulist)[i] <- sampleid
  }
  return(seulist)
}

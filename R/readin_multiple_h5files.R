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

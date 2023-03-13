#' Add mitochondrial information to a Seurat object
#'
#' @param seuratobject single cell data in seurat format
#'
#' @return
#' @export

add_mt <- function(seuratobject) {
  seuratobject[["percent.mt"]] <- Seurat::PercentageFeatureSet(seuratobject, pattern = "^MT-")
  return(seuratobject)
}

#' Add mitochondrial percentage to multiple seurat objects
#'
#' @param seuratlist list of seurat objects
#'
#' @return seurat object with mitochondrial percentages added to it
#' @export

add_multiple_mt <- function(seuratlist){
  seuratlist <- lapply(seuratlist, add_mt)
  return(seuratlist)
}

#' Plot the number of molecules, genes, and mitochondrial percentage per seurat object
#'
#' @param seuratobject from seurat
#' @param min_feature_limit the minimum number of genes you want detected in each cell for QC purposes
#' @param min_mol_limit the minimum number of molecules you want detected in each cell for QC purposes
#' @param plot_title whatever you want the plot to be called
#'
#' @return a plot
#' @export

features_counts_mitoperc_dotplot <- function(seuratobject, min_feature_limit, min_mol_limit, plot_title){
  seuratobject@meta.data %>%
    ggplot2::ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
    ggplot2::geom_point(size = 1) +
    ggplot2::labs(x = "Number of Molecules Detected\nPer Cell", y = "Number of Genes Detected\nPer Cell", color = "MT%") +
    ggplot2::theme_classic() +
    ggplot2::geom_hline(yintercept = min_feature_limit, linetype = 2) +
    ggplot2::geom_vline(xintercept = min_mol_limit, linetype = 2) +
    ggplot2::ggtitle(plot_title)
}

#' Filter one single cell dataset based on molecule count, feature count, and mitochondrial percentage.
#'
#' @param seuratobject the seurat data
#' @param minimum_counts The minimum number of molecules you want per cell
#' @param minimum_features The minimum number of genes you want per cell
#' @param max_mitochondrialperc The maximum mitochondrial percentage
#'
#' @return filtered seurat object
#' @export

qc_filtering <- function(seuratobject,minimum_counts,minimum_features,max_mitochondrialperc){
  subset(seuratobject, nCount_RNA > minimum_counts & nFeature_RNA > minimum_features & percent.mt < max_mitochondrialperc)
  }

#' Filter single cell data based on molecule count, feature count, and mitochondrial percentage.
#'
#' @param seuratlist A list of seurat objects
#' @param minimum_counts The minimum number of molecules you want per cell
#' @param minimum_features The minimum number of genes you want per cell
#' @param max_mitochondrialperc The maximum mitochondrial percentage
#'
#' @return filtered seurat object
#' @export

singlecell_QC_filtering <- function(seuratlist,minimum_counts,minimum_features,max_mitochondrialperc){
  seuratlist_filtered <- lapply(seuratlist, qc_filtering, minimum_counts=minimum_counts,minimum_features=minimum_features,max_mitochondrialperc=max_mitochondrialperc)
  return(seuratlist_filtered)
}

#' Filter a single cell dataset based on molecule count, feature count, and mitochondrial percentage.
#'
#' @param seuratobject A single cell seurat object
#' @param minimum_counts The minimum number of molecules you want per cell
#' @param minimum_features The minimum number of genes you want per cell
#' @param max_mitochondrialperc The maximum mitochondrial percentage
#'
#' @return filtered seurat object
#' @export

onesample_singlecell_QC_filtering <- function(seuratobject,minimum_counts,minimum_features,max_mitochondrialperc){
  seuratobject_filtered <- qc_filtering(seuratobject,minimum_counts=minimum_counts,minimum_features=minimum_features,max_mitochondrialperc=max_mitochondrialperc)
  return(seuratobject_filtered)
}


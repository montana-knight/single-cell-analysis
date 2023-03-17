#' Analyze data for the appropriate number of clusters
#'
#' @param seuratobject your normalized, transformed, etc... seurat object
#' @param variable.selection.method to select your highly variable features
#' @param variable.selection.nfeatures how many variable features you want to pull out
#' @param max.number.clusters.to.check the maximum number of components/clusters you want to look at
#' @param jackstraw.approach true or false to run a jackstraw analysis (takes some time)
#'
#' @return a list with information from both a principal component and jackstraw approach
#' @export

investigate_cluster_number <- function(seuratobject, variable.selection.method,
                                 variable.selection.nfeatures, max.number.clusters.to.check,
                                 jackstraw.approach) {

  # find the most variable features in the dataset to help focus on what is actually
  # different in your dataset and better define clusters
  seuratobject.vf <- FindVariableFeatures(seuratobject, selection.method = "vst", nfeatures = 2000)

  # define all genes
  all.genes <- rownames(seuratobject.vf)

  # scale data to help with extreme counts
  seuratobject.sAll <- ScaleData(seuratobject.vf, features = all.genes)
  seuratobject.sVF <- ScaleData(seuratobject.vf)

  # clustering methods

  # pca cluster
  seuratobject.pAll <- RunPCA(seuratobject.sAll, features = all.genes) # this takes a minute or so...
  seuratobject.pVF <- RunPCA(seuratobject.sVF)

  # cluster examination plots
  pca_dimension_heatmaps <- list()
  for(i in 1:max.number.clusters.to.check){
    pca_dimension_heatmaps[[i]] <- DimHeatmap(seuratobject.pVF, dims = i, cells = 500, balanced = TRUE)
  }

  pca_elbow_plot <- ElbowPlot(seuratobject.pVF)

  pca_clustering_results <- list("pca_varfeat_seurat_object" = seuratobject.pVF,
                                 "dimension_heatmaps" = pca_dimension_heatmaps,
                                 "elbow_plot" = pca_elbow_plot)

  if(jackstraw.approach == TRUE){

    # jackstraw approach

    seuratobject.js <- JackStraw(seuratobject.pVF, num.replicate = 100)
    seuratobject.js <- ScoreJackStraw(seuratobject.js, dims = 1:20)

    jackstrawplot <- JackStrawPlot(seuratobject.js, dims = 1:max.number.clusters.to.check)

    jackstraw_clustering_results <- list("jackstraw_seurat_object" = seuratobject.js,
                                         "jackstraw_plot" <- jackstrawplot)

    clustering_number_investigation_results <- list(pca_clustering_results, jackstraw_clustering_results)

    return(clustering_number_investigation_results)

  } else {

    return(pca_clustering_results)

  }

}

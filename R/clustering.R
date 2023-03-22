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
    pca_dimension_heatmaps[[i]] <- DimHeatmap(seuratobject.pVF, dims = i, cells = 500, balanced = TRUE, fast = FALSE, combine = TRUE)
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
                                         "jackstraw_plot" = jackstrawplot)

    clustering_number_investigation_results <- list("pca_results" = pca_clustering_results, "jackstraw_results" = jackstraw_clustering_results)

    return(clustering_number_investigation_results)

    } else {

      return(pca_clustering_results)

      }
}


#' Apply dimension reduction with Seurat
#'
#' @param seuratobject your normalized, etc... seurat object
#' @param features_to_use either variable or all, for variable features or all genes
#' @param variable.selection.method if using variable
#' @param variable.selection.nfeatures number of features if using variable
#' @param number_of_dimensions number of pca dimensions to use to help estimate clusters
#' @param resolution higher will give more clusters, lower less
#' @param tsne true or false to do a tsne analysis
#' @param umap true or false to do a umap analysis
#'
#' @return
#' @export

sc_clustering <- function(seuratobject, features_to_use, variable.selection.method,
                          variable.selection.nfeatures, number_of_dimensions,
                          resolution, tsne, umap) {
  if(tsne == TRUE || umap == TRUE){

    if(features_to_use == "variable"){

      seuratobject <- FindVariableFeatures(seuratobject, selection.method = variable.selection.method, nfeatures = variable.selection.nfeatures)
      all.genes <- rownames(seuratobject)
      seuratobject <- ScaleData(seuratobject)
      seuratobject <- RunPCA(seuratobject)
      seuratobject <- FindNeighbors(seuratobject, dims = 1:number_of_dimensions)
      seuratobject <- FindClusters(seuratobject, resolution = resolution)

      if(tsne == TRUE && umap == FALSE){
        seuratobject <- RunTSNE(seuratobject, dims = 1:number_of_dimensions)
        tsne_plot <- DimPlot(seuratobject, reduction = "tsne") + coord_equal()
        tsne_analysis <- list("tsne_seurat_object" = seuratobject, "tsne_plot" = tsne_plot)
        return(tsne_analysis)
      } else if(umap == TRUE && tsne == FALSE){
        seuratobject <- RunUMAP(seuratobject, dims = 1:number_of_dimensions)
        umap_plot <- DimPlot(seuratobject, reduction = "umap") + coord_equal()
        umap_analysis <- list("umap_seurat_object" = seuratobject, "umap_plot" = umap_plot)
        return(umap_analysis)
      } else if(umap == TRUE && tsne == TRUE){
        seuratobject <- RunTSNE(seuratobject, dims = 1:number_of_dimensions)
        tsne_plot <- DimPlot(seuratobject, reduction = "tsne") + coord_equal()
        #tsne_analysis <- list("tsne_seurat_object" = seuratobject, "tsne_plot" = tsne_plot)
        seuratobject <- RunUMAP(seuratobject, dims = 1:number_of_dimensions)
        umap_plot <- DimPlot(seuratobject, reduction = "umap") + coord_equal()
        umap_analysis <- list("umap_seurat_object" = seuratobject, "umap_plot"= umap_plot)
        both_analysis <- list("cluster_seurat_object" = seuratobject, "tsne_plot"= tsne_plot, "umap_plot" = umap_plot)
        return(both_analysis)
      }

    } else if(features_to_use == "all"){

      all.genes <- rownames(seuratobject)
      seuratobject <- ScaleData(seuratobject, features = all.genes)
      seuratobject <- RunPCA(seuratobject,features=all.genes)
      seuratobject <- FindNeighbors(seuratobject, dims = 1:number_of_dimensions)
      seuratobject <- FindClusters(seuratobject, resolution = resolution)

      if(tsne == TRUE && umap == FALSE){
        seuratobject <- RunTSNE(seuratobject, dims = 1:number_of_dimensions)
        tsne_plot <- DimPlot(seuratobject, reduction = "tsne") + coord_equal()
        tsne_analysis <- list("tsne_seurat_object" = seuratobject, "tsne_plot" = tsne_plot)
        return(tsne_analysis)
      } else if(umap == TRUE && tsne == FALSE){
        seuratobject <- RunUMAP(seuratobject, dims = 1:number_of_dimensions)
        umap_plot <- DimPlot(seuratobject, reduction = "umap") + coord_equal()
        umap_analysis <- list("umap_seurat_object" = seuratobject, "umap_plot" = umap_plot)
        return(umap_analysis)
      } else {
        seuratobject <- RunTSNE(seuratobject, dims = 1:number_of_dimensions)
        tsne_plot <- DimPlot(seuratobject, reduction = "tsne") + coord_equal()
        #tsne_analysis <- list("tsne_seurat_object" = seuratobject, "tsne_plot"= tsne_plot)
        seuratobject <- RunUMAP(seuratobject, dims = 1:number_of_dimensions)
        umap_plot <- DimPlot(seuratobject, reduction = "umap") + coord_equal()
        umap_analysis <- list("umap_seurat_object" = seuratobject, "umap_plot" = umap_plot)
        both_analysis <- list("cluster_seurat_object" = seuratobject, "tsne_plot"= tsne_plot, "umap_plot" = umap_plot)
        return(both_analysis)
      }

    } else{
      return("Error: invalid feature parameter")
    }

  } else {
    return("Error: Need to choose a reduction technique")
  }

}


#' Apply single cell clustering to an entire list of seurat objects
#'
#' @param seurat_list list of seurat objects
#' @param features_to_use all or variable
#' @param variable.selection.method variable selection method
#' @param variable.selection.nfeatures number of features to use in variable selection
#' @param number_of_dimensions number of dimensions to use for clustering
#' @param resolution clustering resolution
#' @param tsne true or false for tsne reduction
#' @param umap true or false for uman
#'
#' @return
#' @export

clustering_seurat_list <- function(seurat_list, features_to_use, variable.selection.method,
                                   variable.selection.nfeatures, number_of_dimensions,
                                   resolution, tsne, umap){

  clustered_seurat_list <- lapply(seurat_list, sc_clustering, features_to_use, variable.selection.method,
                                  variable.selection.nfeatures, number_of_dimensions,
                                  resolution, tsne, umap)

  return(clustered_seurat_list)

}


#' Annotate the clusters for a seurat object
#'
#' @param seuratobject the object to annotate
#' @param exhaustion_set whether to include the exhaustion signatures
#' @param reduction_type tsne or umap
#'
#' @return plot with clusters defined
#' @export


sctype_annotation <- function(seuratobject, exhaustion_set, reduction_type){

  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"

  if(exhaustion_set == TRUE){
    gs_list <- singlecellanalyzR::sctype_gs_list
  } else{
    gs_list <- singlecellanalyzR::gs_list
  }

  es.max = sctype_score(scRNAseqData = seuratobject[["RNA"]]@scale.data, scaled = TRUE,
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

  cL_resutls = do.call("rbind", lapply(unique(seuratobject@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seuratobject@meta.data[seuratobject@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seuratobject@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  print(sctype_scores[,1:3])

  seuratobject@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,];
    seuratobject@meta.data$customclassif[seuratobject@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }

  annotated_plot <- DimPlot(seuratobject, reduction = reduction_type, label = TRUE, repel = TRUE, group.by = 'customclassif')


  # cL_resutls=cL_resutls[order(cL_resutls$cluster),]; edges = cL_resutls; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL
  # #
  # # # prepare nodes
  # nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c();
  # ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
  # for (i in 1:length(unique(cL_resutls$cluster))){
  #    dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
  #  }
  #  nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
  #  files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
  #  nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]
  # #
  #  mygraph <- graph_from_data_frame(edges, vertices=nodes)
  # #
  # # # Make the graph
  #  gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) +
  #    geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  #    theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")
  # #
  #  multiplot_annot <- scater::multiplot(DimPlot(seuratobject, reduction = reduction_type, label = TRUE, repel = TRUE, cols = ccolss), gggr, cols = 2)
  # #
   output <- list(seuratobject, annotated_plot)
   return(output)
}

#' Apply cluster annotion with sc type on all objects in a list
#'
#' @param seuratlist a list of seurat objects
#' @param exhaustion_set true or false to use the exhaustion data
#' @param reduction_type tsne or umap
#'
#' @return a list with the new seurat objects, and plots
#' @export

sctype_list_annotation <- function(seuratlist, exhaustion_set, reduction_type){
  new_seurat_list <- lapply(seuratlist, sctype_annotation, exhaustion_set, reduction_type)
  return(new_seurat_list)
}

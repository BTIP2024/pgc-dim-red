# functions for the pre-processing of scRNAseq

# creating the seurat object from .h5 and .gz files
write_seurat_object <- function(input) {
   if(grepl('.gz', input)){
      utils::untar(input, files = NULL, list = FALSE, exdir = "sample")
      new <- list.files(path = "./sample/", pattern = "*.gz", recursive = TRUE, full.names = T) #WORKED
      
      file.copy(from = new, to = ".", overwrite = T)
      
      obj <- Seurat::Read10X(data.dir = ".")
   } else if(grepl('.loom', input)) {
      obj <- SeuratDisk::Connect(filename = input, mode = "r")
      obj <- Seurat::as.Seurat(obj)
   } else if(grepl('.h5ad', input)){
      SeuratDisk::Convert(input, dest = "h5seurat", overwrite = TRUE)
      h5seurat <- list.files(path = ".", pattern = ".h5seurat", recursive = TRUE, full.names = T)
      obj <- SeuratDisk::LoadH5Seurat(h5seurat)
   } else {
      h5_data <- Seurat::Read10X_h5(filename = input, use.names = TRUE, unique.features = TRUE)
      
      if(class(h5_data) == "list") {
         cts <- h5_data$'Gene Expression'
         obj <- Seurat::CreateSeuratObject(counts = cts, project = "Seurat", min.cells = 3, min.features = 200)
      } else {
         obj <- Seurat::CreateSeuratObject(counts = h5_data, project = "Seurat", min.cells = 3, min.features = 200)
      } 
   }
   saveRDS(obj, file = "seurat_object.rds")
}


# quality control

qc_seurat <- function(input) {
   seurat_obj <- readRDS(input)
   seurat_obj <- Seurat::UpdateSeuratObject(seurat_obj)
   seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^MT-")
   
   vplot <- Seurat::VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + ggplot2::theme(legend.position = "none") 
   Seurat::FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
   
   ggplot2::ggsave(vplot, file = "violin_plot.png", width = 15)
   
   vplot1 <- Seurat::VlnPlot(seurat_obj, features = "nFeature_RNA") + ggplot2::theme(legend.position = "none")
   vplot2 <- Seurat::VlnPlot(seurat_obj, features = "nCount_RNA") + ggplot2::theme(legend.position = "none")
   vplot3 <- Seurat::VlnPlot(seurat_obj, features = "percent.mt") + ggplot2::theme(legend.position = "none")
   
   library(plotly)
   
   ggvplot1 <- ggplotly(vplot1)
   ggvplot2 <- ggplotly(vplot2)
   ggvplot3 <- ggplotly(vplot3)
   
   annotations = list(list(x = 0.15, y = 1, text = "nFeature_RNA", xref = "paper", yref = "paper", xanchor = "center", yanchor = "bottom", showarrow = FALSE), list(x = 0.5, y = 1, text = "nCount_RNA", xref = "paper", yref = "paper", xanchor = "center", yanchor = "bottom", showarrow = FALSE), list(x = 0.85, y = 1, text = "percent.mt", xref = "paper", yref = "paper", xanchor = "center", yanchor = "bottom", showarrow = FALSE))
   
   threeplots <- plotly::subplot(ggvplot1, ggvplot2, ggvplot3)
   threeplots %>% plotly::layout(title = "Violin Plots", annotations = annotations)
   
   htmltools::save_html(threeplots, file = "violin_plots.html")
   
   splot1 <- Seurat::FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") + ggplot2::theme(legend.position = "none")
   splot2 <- Seurat::FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggplot2::theme(legend.position = "none")

   ggsplot1 <- ggplotly(splot1)
   ggsplot2 <- ggplotly(splot2)
   
   annotations = list(list(x = 0.15, y = 1, text = "-0.05", xref = "paper", yref = "paper", xanchor = "center", yanchor = "bottom", showarrow = FALSE), list(x = 0.5, y = 1, text = "0.95", xref = "paper", yref = "paper", xanchor = "center", yanchor = "bottom", showarrow = FALSE))
   
   scatterplots <- plotly::subplot(ggsplot1, ggsplot2)
   scatterplots %>% plotly::layout(title = "Scatter Plots", annotations = annotations)
   
   htmltools::save_html(scatterplots, file = "scatter_plots.html")
   
   seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
   saveRDS(seurat_obj, file = "afterQC_seurat.rds")
}


# data normalization 

normalization_seurat <- function(input) {
   seurat_obj <- readRDS(input)
   seurat_obj <- Seurat::NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
   saveRDS(seurat_obj, file = "afternormalization.rds")
}


# identifying variable features

variable_features <- function(input){
   normalized <- readRDS(input)
   normalized <- Seurat::FindVariableFeatures(normalized, selection.method = "vst", nfeatures = 2000)
   
   top10 <- head(Seurat::VariableFeatures(normalized), 10)
   
   plot1 <- Seurat::VariableFeaturePlot(normalized)
   plot2 <- Seurat::LabelPoints(plot = plot1, points = top10, repel = TRUE)
   finalplots <- plot1 + plot2
   
   ggplot2::ggsave(finalplots, filename = "variable_features_plot.png", width = 15, height = 8, dpi = 300)
   saveRDS(top10, file = "top10_variablefeatures.rds")
   saveRDS(normalized, file = "variablefeatures.rds")
}

# scaling
scaling_seurat <- function(input){
   for_scaling <- readRDS(input)
   all.genes <- rownames(for_scaling)
   for_scaling <- Seurat::ScaleData(for_scaling, features = all.genes)
   saveRDS(for_scaling, file = "after_scaling.rds")
}

# pca
pca_seurat <- function(input){
   for_pca <- readRDS(input)
   for_pca <- Seurat::RunPCA(for_pca, features = Seurat::VariableFeatures(object = for_pca))
   print(for_pca[["pca"]], dims = 1:5, nfeatures = 5)
   
   image <- Seurat::VizDimLoadings(for_pca, dims = 1:2, reduction = "pca")
   ggplot2::ggsave(image, file = "pca_results.png", width = 15, height = 10)
   
   image1 <- Seurat::DimPlot(for_pca, reduction = "pca", label = TRUE)
   ggplot2::ggsave(image1, file = "pca_plot_unlabeled.png", width = 15, height = 10)
   
   image2 <- Seurat::DimHeatmap(for_pca, dims = 1, cells = 500, balanced = TRUE)
   ggplot2::ggsave(image2, file = "heatmap.png", width = 12, height = 12)
   
   image3 <- Seurat::DimHeatmap(for_pca, dims = 1:10, cells = 500, balanced = TRUE)
   ggplot2::ggsave(image3, file = "heatmap_multiple.png", width = 20, height = 20)
   
   image4 <- Seurat::ElbowPlot(for_pca)
   ggplot2::ggsave(image4, file = "elbowplot.png", width = 10)
}

# non-linear dimensionality reduction
clusters_seurat <- function(input){
   clustering <- readRDS(input)
   clustering <- Seurat::RunPCA(clustering, features = Seurat::VariableFeatures(object = clustering))
   image4 <- Seurat::FindNeighbors(clustering, dims= 1:15)
   image4 <- Seurat::FindClusters(image4, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
   saveRDS(image4, file = "clustersperresolution.rds")
   
   image4 <- Seurat::DimPlot(image4, group.by = "RNA_snn_res.0.3", label = TRUE)
   
   ggplot2::ggsave(image4, file = "dimplot.png", width = 12, height = 10)
}

tsne_seurat <- function(input){
   for_tsne <- readRDS(input)
   for_3d <- for_tsne
   
   tsne <- Seurat::RunTSNE(for_tsne, dims = 1:10, dim.embed =2, label = TRUE)
   tsne <- Seurat::DimPlot(tsne, reduction = "tsne")
   ggplot2::ggsave(tsne, file = "tsne_seurat.png", width = 10, height = 10)
   
   
   #for 3D plots
   tsne_3d <- Seurat::RunTSNE(for_3d, dims = 1:10, dim.embed = 3)
   
   tsne_1 <- tsne_3d[["tsne"]]@cell.embeddings[,1]
   tsne_2 <- tsne_3d[["tsne"]]@cell.embeddings[,2]
   tsne_3 <- tsne_3d[["tsne"]]@cell.embeddings[,3]
   
   plot.data <- Seurat::FetchData(object = tsne_3d, vars = c("tSNE_1", "tSNE_2", "tSNE_3", "seurat_clusters"))
   
   plot.data$label <- paste(rownames(plot.data))
   
   plotin3d <- plotly::plot_ly(data = plot.data, 
           x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
           color = ~seurat_clusters,
           type = "scatter3d", 
           mode = "markers", 
           marker = list(size = 5, width=2), # controls size of points
           text=~label, #This is that extra column we made earlier for which we will use
           hoverinfo="text")
   
   htmltools::save_html(plotin3d, file = "tsne_3dplot.html")

}

umap_seurat <- function(input){
   for_umap <- readRDS(input)
   
   umap <- Seurat::RunUMAP(for_umap, dims = 1:10, n.components = 3L)
   umap <- Seurat::DimPlot(umap, reduction = "umap")
   ggplot2::ggsave(umap, file = "tsne_seurat.png", width = 10, height = 10)
   
   plot.data <- Seurat::FetchData(object = for_umap, vars = c("umap_1", "umap_2", "umap_3", "seurat_clusters"))
   
   plot.data$label <- paste(rownames(plot.data))
   
   fig <- plot_ly(data = plot.data, 
                  x = ~umap_1, y = ~umap_2, z = ~umap_3, 
                  color = ~seurat_clusters,
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 5, width=2), # controls size of points
                  text=~label, 
                  hoverinfo="text")
   
   htmltools::save_html(fig, file = "umap_3dplot.html")
}



# ADDITIONAL ONLY

gene_expression <- function(input, gene){
   for_tsne <- readRDS(input)
   gene_input <- gene
   tsne_3d <- Seurat::RunTSNE(for_tsne, dims = 1:10, dim.embed = 3)
   
   tsne_1 <- tsne_3d[["tsne"]]@cell.embeddings[,1]
   tsne_2 <- tsne_3d[["tsne"]]@cell.embeddings[,2]
   tsne_3 <- tsne_3d[["tsne"]]@cell.embeddings[,3]
   
   plotting.data <- Seurat::FetchData(object = pbmc, vars = c("tSNE_1", "tSNE_2", "tSNE_3", gene_input))
   
   plotting.data$changed <- ifelse(test = plotting.data$gene_input <1, yes = plotting.data$gene_input, no = 1)
   
   plotting.data$label <- paste(rownames(plotting.data)," - ", plotting.data$gene_input, sep="")
   
   ggplot2::plot_ly(data = plotting.data, 
                    x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
                    color = ~changed,
                    opacity = .5,
                    colors = c('darkgreen', 'red'), 
                    type = "scatter3d", 
                    mode = "markers",
                    marker = list(size = 5, width=2), 
                    text=~label,
                    hoverinfo="text"
   )
}
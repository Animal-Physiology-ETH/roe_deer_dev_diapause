#' @title anno_legend
#' @description
#' graphic param for annotation legend
#' @import ComplexHeatmap
#'

anno_legend = function(fonts_min){
  return(list(title_gp = gpar(fontsize = fonts_min+2,
                              fontface = "bold"),
              labels_gp = gpar(fontsize = fonts_min)))
}

#' @title embryo_hm_annot
#' @description
#' heatmap annotations for embryos
#' @import ComplexHeatmap
#'
#'
#' @param expr_mat numeric matrix containing expression data
#' @param clusterEmbryos data frame containing column with embryo identifier (embryo_ID_column) and cluster assignment
#' @param ClusterIDcolumnE column name in clusterEmbryos with names for embryo / sample clusters
#' @param embryo_annot data frame containing annotations for embryo data
#' @param embryo_ID_column column name with embryo / sample ID variable (must be identical to either column or row names of expr_mat)
#' @param cell_col_name column name in embryo_annot containing log2-transformed cell numbers
#' @param morph_col_name column name in embryo_annot containing log2-transformed cell numbers
#' @param colors_cells color vector with 3 elements for cell numbers
#' @param colors_morph data frame containing color assignments for morphologies
#' @param colors_cluster data frame containing color assignments for cluster
#' @param Embryoskm number of embryo / sample clusters for kmeans clustering
#' @param Embryoskm_repeats number of repeats for kmeans clustering for embryo / sample clusters
#' @param RC direction of annotation ("col" or "row")
#' @param show_lg display annotation legend (T/F)
#' @param anno_size size of annotation (unit())
#' @param base_size base size for fonts
#' @param return_cellno return annotation for cell number? (T/F)
#' @param return_morph return annotation for morphology? (T/F)
#' @param return_cluster return annotation for clustering? (T/F)
#' @export

embryo_hm_annot = function(expr_mat,
                           clusterEmbryos = clustering$Embryos,
                           ClusterIDcolumnE = "ClusterNo",
                           embryo_annot = embryo_info,
                           embryo_ID_column = "E_ID",
                           cell_col_name = "Log2Cells",
                           morph_col_name = "Morphology",
                           colors_cells = (colorscheme %>% filter(Class == "Cells"))$Color,
                           colors_morph = colorscheme %>% filter(Class == "Morphology"),
                           colors_cluster = (colorscheme %>% filter(Class == "EmbryoCluster")),
                           RC = "col",
                           show_lg = TRUE,
                           anno_size = unit(2, "mm"),
                           show_name = FALSE,
                           base_size = 8,
                           return_cellno = F,
                           return_morph = F,
                           return_cluster = F){

  df = embryo_annot[match(colnames(expr_mat), embryo_annot[[embryo_ID_column]]),]
  return_annos = c()

  cellno  = HeatmapAnnotation(Log2Cells= df[[cell_col_name]],
                                   annotation_label = "Cell #\n(log2)",
                                   which = RC,
                                   annotation_legend_param = anno_legend(fonts_min = base_size),
                                   col =  list(Log2Cells= circlize::colorRamp2(c(min(df[[cell_col_name]]),
                                                                 median(df[[cell_col_name]]),
                                                               max(df[[cell_col_name]])),
                                                               colors_cells)),
                                   simple_anno_size = anno_size,
                                   show_annotation_name = show_name)
  if (return_cellno == T) {
    return(cellno)
  }

  morph = HeatmapAnnotation(Morphology=df[[morph_col_name]],
                                 annotation_label = "Morphology",
                                 which = RC,
                                 annotation_legend_param = anno_legend(fonts_min = base_size),
                                 col =  list(Morphology = setNames(colors_morph$Color, colors_morph$Item)),
                                 simple_anno_size = anno_size,
                                 show_annotation_name = show_name)
  if (return_morph == T) {
    return(morph)
  }
  clusterEmbryos
  E_clustering = clusterEmbryos[match(colnames(expr_mat), clusterEmbryos[[embryo_ID_column]]),]
  cluster = HeatmapAnnotation(EmbryoclusterColors= setNames(E_clustering[[ClusterIDcolumnE]],E_clustering[[embryo_ID_column]]),
                                   which = RC,
                                   #annotation_legend_param = anno_legend(fonts_min = fsize),
                                   show_legend = FALSE,
                                   col = list(EmbryoclusterColors = setNames(colors_cluster$Color, colors_cluster$Item)),
                                   simple_anno_size = anno_size,
                                   show_annotation_name = FALSE)
  if (return_cluster == T) {
    return(cluster)
  }
}

#' @title genes_hm_annot
#' @description
#' heatmap annotations for embryos
#' @import ComplexHeatmap
#' @param expr_mat numeric matrix containing expression data
#' @param clusterGenes data frame containing column with Gene identifier (GeneName) and cluster assignment (ABC)
#' @param ClusterIDcolumnG column name in clusterGenes with names for gene clusters
#' @param gene_ID_column column name with gene ID (must be identical to either column or row names of expr_mat).
#' @param gene_annot data frame containing annotations for genes (also defines which genes to display)
#' @param grouping_var variable in gene_annot containing the groups to be displayed (names must be found in coloring)
#' @param sorting_var variable in gene_annot for which legend should be sorted
#' @param colors data frame with colorscheme
#' @param RC direction of annotation ("col" or "row")
#' @param show_lg display annotation legend (T/F)
#' @param anno_size size of annotation (unit())
#' @param show_name show annotation name (T/F)
#' @param base_size base size for fonts


genes_hm_annot = function(expr_mat,
                          clusterGenes = clustering$Genes,
                          ClusterIDcolumnG = "ClusterNo",
                          gene_ID_column = "GeneName",
                          gene_annot = markergenes,
                          grouping_var = "Compartment",
                          sorting_var = "Compartment_No",
                          colors = colorscheme,
                          RC = "row",
                          show_lg = TRUE,
                          anno_size = unit(4, "mm"),
                          show_name = FALSE,
                          base_size = 8){
  gene_annot = gene_annot[(gene_annot[[gene_ID_column]] %in% rownames(expr_mat)),]
  gene_annot = gene_annot %>% select(c(gene_ID_column, grouping_var, sorting_var))
  gene_annot = gene_annot[order(gene_annot[[sorting_var]]),]
  gene_annot[[grouping_var]] <- factor(gene_annot[[grouping_var]], levels=as.factor(unique(gene_annot[[grouping_var]])))

  colors_ = colors %>% filter(Class == grouping_var)

  anno_group  = HeatmapAnnotation(MarkerGroup= gene_annot[[grouping_var]],
                                  annotation_label = grouping_var,
                                  which = RC,
                                  annotation_legend_param = list(title_gp = gpar(fontsize = base_size+2,
                                                                                 fontface = "bold"),
                                                                 labels_gp = gpar(fontsize = base_size),
                                                                 MarkerGroup = list(at  = levels(gene_annot[[grouping_var]]))),
                                  col = list(MarkerGroup = setNames(colors_$Color, colors_$Item)),
                                  simple_anno_size = anno_size,
                                  show_annotation_name = show_name)

  return(c(anno_group))
}

#stack heatmap vertically %v%

#' @title figure2
#' @description
#' plots heatmap using ComplexHeatmap
#'
#' @param expr_mat numeric matrix containing expression data
#' @param gene_as_rows are genes as rows? (T/F) --> if false matrix will by transposed
#' @param compute_clusters_new re-compute clusters for embryos / samples and genes
#' @param seed seed used for kmeans clustering
#'
#' @param embryo_annot data frame containing annotations for embryo data
#' @param embryo_ID_column column name with embryo / sample ID variable (must be identical to either column or row names of expr_mat)
#' @param clusterEmbryos data frame with cluster assignment for embryos / samples (ignored if compute_clusters_new = T)
#' @param ClusterIDcolumnE column name in clusterEmbryos with names for embryo / sample clusters
#' @param Embryoskm number of embryo / sample clusters for kmeans clustering
#' @param Embryoskm_repeats number of repeats for kmeans clustering for embryo / sample clusters
#'
#' @param gene_annot data frame containing annotations for genes (also defines which genes to display)
#' @param gene_ID_column column name with gene ID (must be identical to either column or row names of expr_mat).
#' Order of elements in gene_ID_column will be maintained inside each gene cluster.
#' @param clusterGenes data frame with cluster assignment for genes (ignored if compute_clusters_new = T)
#' @param ClusterIDcolumnG column name in clusterGenes with names for gene clusters
#' @param Genekm number of genee clusters for kmeans clustering
#' @param Genekm_repeats number of repeats for kmeans clustering for gene clusters
#' @param grouping_var variable in gene_annot containing the groups to be displayed (names must be found in coloring)
#'
#' @param base_size base size for text (pt)
#' @param colors dataframe with color assignement
#' @param scaling was data scaled minmax or zscore (will change colors of heatmap)
#' @param show_expr_lvls should absolute expression levels be shown (T/F)? If T, provide unscaled_mat
#' @param value_name name for heatmap variable
#'
#' @param filter_data should data be filtered (T/F)? If T, provide unscaled_mat, min_TMM, and min_s.
#' @param unscaled_mat dataframe with unscaled expression values (identical format as expr_mat)
#' @param min_TMM min TMM value that needs to be observed in at least min_s samples for a gene to be shown
#' @param min_s min number of samples that need to display at least TMM value for a gene to be shown
#'
#' @import ComplexHeatmap
#' @import dplyr
#' @import viridis
#' @importFrom circlize colorRamp2
#' @export

figure2 = function(expr_mat = exprmat_TMM_minmax,
                          gene_as_rows = T,
                          compute_clusters_new = F,
                          seed = 123,

                          embryo_annot = embryo_info,
                          embryo_ID_column = "E_ID",
                          clusterEmbryos = clustering$Embryos,
                          ClusterIDcolumnE = "ClusterNo",
                          Embryoskm = 4,
                          Embryoskm_repeats = 800,

                          gene_annot = markergenes,
                          gene_ID_column = "GeneName",
                          clusterGenes = clustering$Genes,
                          ClusterIDcolumnG = "ClusterNo",
                          Genekm = 6,
                          Genekm_repeats = 800,
                          grouping_var = "Compartment",

                          base_size = 8,
                          colors = colorscheme,
                          scaling = "minmax",
                          show_expr_lvls = T,
                          value_name = "scaled\nTMM\ncounts",

                          filter_data = F,
                          unscaled_mat  = exprmat_TMM,
                          min_TMM = 1,
                          min_s = 3){
  if (gene_as_rows != T) {
    expr_mat = t(expr_mat)
  }

  #make sure all samples are listed in embryo info
  embryo_annot = embryo_annot[embryo_annot[[embryo_ID_column]] %in% colnames(expr_mat),]
  expr_mat = expr_mat[,colnames(expr_mat) %in% embryo_annot[[embryo_ID_column]]]

  gene_annot = gene_annot[gene_annot[[gene_ID_column]] %in% rownames(expr_mat),]
  expr_mat = expr_mat[rownames(expr_mat) %in% gene_annot[[gene_ID_column]],]

  #get order from gene_annot and match expr_mat to that order
  expr_mat = expr_mat[match(gene_annot[[gene_ID_column]], rownames(expr_mat)),]

  # match unscaled data
  unscaled_matching = unscaled_mat[rownames(unscaled_mat) %in% rownames(expr_mat),]
  unscaled_matching = unscaled_matching[match(gene_annot[[gene_ID_column]], rownames(unscaled_matching)),]

  if (filter_data == T) {
    filter_selected = unscaled_matching
    if (all(names(filter_selected == expr_mat))) {
        expr_mat_f =  expr_mat[rowSums(
        as.data.frame(filter_selected > min_TMM) %>%
          mutate(across(everything(), ~+as.logical(.x)))) > min_s,]
        expr_mat = expr_mat_f
    }else{
      stop("Please provide a filter df with same dimensions and names as expression df!")
    }
  }


  if (show_expr_lvls == T) {
    unscaled_matching = log2(unscaled_matching+0.1)
    expr_lvls = HeatmapAnnotation(which = "row",
                                  name = "unscaled TMM",
                                  annotation_name_gp = grid::gpar(fontsize = base_size),
                                  show_legend = F,
                                  annotation_label = "log2(TMM)",
                                  TMM = anno_boxplot(unscaled_matching, height = unit(4, "cm")),
                                  border = F)
  } else{
    expr_lvls = NULL
  }

  if (compute_clusters_new == T) {
    cluster_new = clustering_kmeans(exprmat = expr_mat,
                      selection = rownames(expr_mat),
                      embryo_i = embryo_annot,
                      seed = seed,
                      rowkm = Genekm,
                      columnkm = Embryoskm,
                      rowkm_rep = Genekm_repeats,
                      columnkm_rep = Embryoskm_repeats,
                      save =F)
    names(cluster_new$Embryos) = c(embryo_ID_column, ClusterIDcolumnE)
    clusterEmbryos = cluster_new$Embryos
    names(cluster_new$Genes) = c(gene_ID_column, ClusterIDcolumnG)
    clusterGenes = cluster_new$Genes
  }

  embryo_km_title = list(unique(clusterEmbryos[[ClusterIDcolumnE]]))[[1]]
  gene_km_title = list(unique(clusterGenes[[ClusterIDcolumnG]]))[[1]]


  if (scaling == "minmax") {
    hm_cols = colors %>% filter(Class=="HeatmapMM")
    } else{
    hm_cols = colors %>% filter(Class=="HeatmapZ")
  }
  color_fn = circlize::colorRamp2(as.numeric(hm_cols$Item), hm_cols$Color)

  set.seed(seed)
  ht = Heatmap(expr_mat,
               col = color_fn,
               name = value_name,
               row_title = gene_km_title,
               row_km = Genekm,
               row_km_repeats = Genekm_repeats,
               show_row_dend = FALSE,
               cluster_rows = FALSE,
               cluster_row_slices = FALSE,

               column_km = Embryoskm,
               column_title = embryo_km_title,
               column_km_repeats = Embryoskm_repeats,
               show_column_dend = FALSE,
               show_column_names = FALSE,

               row_title_gp = grid::gpar(fontsize = base_size+2, fontface = "bold"),
               row_names_gp =  grid::gpar(fontsize = base_size),
               column_title_gp = grid::gpar(fontsize = base_size+2, fontface = "bold"),

               # heatmap_legend_param = anno_legend(base_size),
               bottom_annotation  = c(embryo_hm_annot(expr_mat = expr_mat,
                                                      embryo_annot = embryo_annot,
                                                      clusterEmbryos = clusterEmbryos,
                                                      ClusterIDcolumnE = ClusterIDcolumnE,
                                                      embryo_ID_column = embryo_ID_column,
                                                      colors_cells = (colors %>% filter(Class == "Cells"))$Color,
                                                      base_size = base_size,
                                                      return_cellno = T),
                                      embryo_hm_annot(expr_mat = expr_mat,
                                                      embryo_annot = embryo_annot,
                                                      clusterEmbryos = clusterEmbryos,
                                                      ClusterIDcolumnE = ClusterIDcolumnE,
                                                      embryo_ID_column = embryo_ID_column,
                                                      colors_morph = colors %>% filter(Class == "Morphology"),
                                                      base_size = base_size,
                                                      return_morph =  T)),
               right_annotation = genes_hm_annot(expr_mat = expr_mat,
                                                 clusterGenes = clusterGenes,
                                                 ClusterIDcolumnG = ClusterIDcolumnG,
                                                 gene_ID_column = gene_ID_column,
                                                 gene_annot = gene_annot,
                                                 grouping_var = grouping_var,
                                                 colors = colors,
                                                 base_size = base_size),
               top_annotation = embryo_hm_annot(expr_mat = expr_mat,
                                                embryo_annot = embryo_annot,
                                                clusterEmbryos = clusterEmbryos,
                                                ClusterIDcolumnE = ClusterIDcolumnE,
                                                embryo_ID_column = embryo_ID_column,
                                                colors_cluster = (colors %>% filter(Class == "EmbryoCluster")),
                                                base_size = base_size,
                                                return_cluster = T),
               left_annotation = expr_lvls,

               show_heatmap_legend = TRUE)
  return(ht)
}










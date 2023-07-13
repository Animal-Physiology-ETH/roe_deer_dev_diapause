#PATHWAWY LIST BASED ON : https://www.nature.com/articles/s41422-021-00592-9

#  Create a heatmap for pluripotency hub-genes
#' @title metadata_rda
#' @description creates an rda objects of all metadata
#' @param pp_pws  list with gene names and pathway association
#' @param gene_ID_column column name with gene ID (must be identical to either column or row names of expr_mat).
#' @param grouping_var variable in pp_pws containing the groups to be displayed
#' @param expr_mat numeric matrix containing expression data
#' @param embryo_annot data frame containing annotations for embryo data
#' @param embryo_ID_column column name with embryo / sample ID variable (must be identical to either column or row names of expr_mat)
#' @param clusterEmbryos data frame with cluster assignment for embryos / samples (ignored if compute_clusters_new = T)
#' @param ClusterIDcolumnE column name in clusterEmbryos with names for embryo / sample clusters
#' @param orderE order of embryos (cluster (precomputed based on figure 2 heatmap) or cells)
#' @param colors dataframe with color assignement
#' @param filter_data should data be filtered (T/F)? If T, provide unscaled_mat, min_TMM, and min_s.
#' @param unscaled_mat dataframe with unscaled expression values (identical format as expr_mat)
#' @param min_TMM min TMM value that needs to be observed in at least min_s samples for a gene to be shown
#' @param min_s min number of samples that need to display at least TMM value for a gene to be shown
#' @param base_size base size for text (pt)
#'
#'
#' @returns heatmap object
#' @import ComplexHeatmap
#' @import dplyr
#' @import viridis
#' @importFrom circlize colorRamp2
#' @export
#'

pppw_heatmap = function(pp_pws = pppw_hubgenes,
                        gene_ID_column = "GeneName",
                        grouping_var = "Pathway",
                        expr_mat = exprmat_TMM_minmax,
                        embryo_annot = embryo_info,
                        embryo_ID_column = "E_ID",
                        clusterEmbryos = clustering$Embryos,
                        ClusterIDcolumnE = "ClusterNo",
                        orderE = c("cluster", "cells"),
                        colors = colorscheme,
                        filter_data = F,
                        unscaled_mat  = exprmat_TMM,
                        min_TMM = 1,
                        min_s = 3,
                        base_size= 8){
  expr_df = as.data.frame(t(expr_mat))

  expr_selected = expr_df %>% select(pp_pws[[gene_ID_column]][pp_pws[[gene_ID_column]] %in% colnames(expr_df)])

  if (filter_data == T) {
    unscaled_matching = unscaled_mat[rownames(unscaled_mat) %in% rownames(expr_mat),]
    unscaled_matching = unscaled_matching[match(pp_pws[[gene_ID_column]], rownames(unscaled_matching)),]
    expr_mat_sel = expr_mat[match(pp_pws[[gene_ID_column]], rownames(expr_mat)),]
    if (all(names(unscaled_matching == expr_mat_sel))) {
      expr_mat_f =  expr_mat_sel[rowSums(
        as.data.frame(unscaled_matching > min_TMM) %>%
          mutate(across(everything(), ~+as.logical(.x)))) > min_s,]
      expr_mat = expr_mat_f
    }else{
      stop("Please provide a filter df with same dimensions and names as expression df!")
    }
  }

  pp_pws = pp_pws %>% filter(GeneName %in% colnames(expr_selected))
  pp_pws = pp_pws %>% filter(!duplicated(pp_pws$GeneName))
  pp_pws$Group <- factor(pp_pws[[grouping_var]], levels=as.factor(unique(pp_pws[[grouping_var]])))

  expr_matching = t(expr_selected)

  embryo_anno = embryo_hm_annot(expr_mat = expr_matching, embryo_annot = embryo_annot, base_size = base_size)

  top_anno = embryo_hm_annot(expr_mat = expr_matching, embryo_annot = embryo_annot, base_size = base_size, return_cluster = TRUE)

  clusterE=clusterEmbryos[match(colnames(expr_matching), clusterEmbryos[[embryo_ID_column]]),]


  if (orderE == "cluster") {
    colsplit = clusterE$ClusterNo
    column_order = order(clusterE$ClusterNo)
    colreorder = F}
  if (orderE == "cells") {
    colsplit = NULL
    column_order = order(embryo_annot$Cells)
    colreorder = F
  }

  hm_cols = colors %>% filter(Class=="HeatmapMM")
  color_fn = circlize::colorRamp2(as.numeric(hm_cols$Item), hm_cols$Color)

  unscaled_data = unscaled_mat[rownames(expr_matching),]
  unscaled_data = log2(unscaled_data+0.1)

  ha = HeatmapAnnotation(which = "row",
                         name = "unscaled TMM",
                         annotation_name_gp = grid::gpar(fontsize = base_size),
                         show_legend = F,
                         annotation_label = "log2(TMM)",
                         border = F,
                         TMM = anno_boxplot(unscaled_data, height = unit(4, "cm")))

  ht = ComplexHeatmap::Heatmap(as.matrix(expr_matching),
                               col = color_fn,
                               row_title_gp = grid::gpar(fontsize = base_size, fontface = "bold"),
                               column_title_gp = gpar(fontsize = base_size, fontface = "bold"),
                               show_column_names = F,
                               row_names_gp =  grid::gpar(fontsize = base_size),
                               column_names_gp =  grid::gpar(fontsize = base_size),

                               name = "scaled\nTMM\ncounts",
                               cluster_row_slices =F,
                               row_split  = pp_pws$Group,
                               row_dend_reorder = F,
                               column_split = colsplit,
                               column_order = column_order,
                               column_dend_reorder = colreorder,
                               bottom_annotation  = embryo_anno,
                               left_annotation = ha,
                               top_annotation = top_anno,
                               show_column_dend = FALSE,
                               show_row_dend = FALSE
  )
  return(ht)
  }


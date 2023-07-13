#' @title plot_pca
#' @param expr_mat numeric matrix containing expression data, colnames (or rownames) should contain values found in embryo_ID_column
#' @param gene_as_rows are genes as rows? (T/F) --> if false matrix will by transposed
#'
#' @param embryo_annot data frame containing annotations for embryo data
#' @param embryo_ID_column column name with embryo / sample ID variable (must be identical to either column or row names of expr_mat)
#'
#' @param clusterEmbryos data frame with cluster assignment for embryos / samples (ignored if compute_clusters_new = T)
#' @param ClusterIDcolumnE column name in clusterEmbryos with names for embryo / sample clusters
#'
#' @param gene_annot data frame containing annotations for genes (also defines which genes to display)
#' @param gene_ID_column column name with gene ID (must be identical to either column or row names of expr_mat).
#' Order of elements in gene_ID_column will be maintained inside each gene cluster.
#' @param grouping_var variable in gene_annot containing the groups to be displayed
#'
#' @param base_size base size for text (pt)
#' @param colorsE vector with colors for embryo clusters
#' @param colorsG vector with colors for genes
#' @param output which plot(s) to return c("ind", biplot","both")
#'
#' @import factoextra
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @importFrom cowplot get_legend

plot_pca = function(expr_mat = exprmat_TMM_minmax,
                    gene_as_rows = T,

                    embryo_annot = embryo_info,
                    embryo_ID_column = "E_ID",

                    clusterEmbryos = clustering$Embryos,
                    ClusterIDcolumnE = "ClusterNo",

                    gene_annot = markergenes,
                    gene_ID_column = "GeneName",
                    grouping_var = "Compartment",

                    colorsE = colorscheme[colorscheme$Class == "EmbryoCluster",][["Color"]],
                    colorsG = colorscheme[colorscheme$Class ==  "Compartment",][["Color"]],
                    base_size = 8,
                    colorlab = "Cluster\nNumber",
                    output = "ind"
                    ){
  if (gene_as_rows != T) {
    expr_mat = t(expr_mat)
  }
  expr_mat = expr_mat[rownames(expr_mat) %in% gene_annot[[gene_ID_column]],]
  expr_mat = expr_mat[,order(colnames(expr_mat))]
  #create same order as embryo info df
  embryo_annot = embryo_annot[order(embryo_annot[[embryo_ID_column]]),]
  #select only embryos that are also annotated
  expr_mat= expr_mat[(colnames(expr_mat) %in% embryo_annot[[embryo_ID_column]]),]

  clusterEmbryos = clusterEmbryos[match(colnames(expr_mat), clusterEmbryos[[embryo_ID_column]]),]

  res.pca = prcomp(t(expr_mat), scale = FALSE)

  pca_ind = factoextra::fviz_pca_ind(res.pca,
                                     addEllipses = TRUE, # Concentration ellipses,
                                     col.ind = clusterEmbryos[[ClusterIDcolumnE]],
                                     pointsize = 1 ,
                                     labelsize = 1,
                                     repel = TRUE,
                                     pointshape = 19,
                                     geom = "point")+
    scale_fill_manual(values = colorsE)+
    scale_color_manual(values = colorsE)+
    labs(color=colorlab)+
    guides(fill = "none",
           color = guide_legend(override.aes = list(shape = 16, label = "", size = 1),
                                direction = "vertical",
                                title.position = "top",
                                title.theme = element_text(size = base_size, face = "bold")))+
    theme_general()+
    theme(legend.position = c(0.05, 0.95),
          legend.direction = "vertical")

  if (output == "ind") {
    return(pca_ind)
  }

  # biplot
  gene_annot = gene_annot[gene_annot[[gene_ID_column]] %in% rownames(expr_mat),]
  gene_annot = gene_annot[match(rownames(expr_mat), gene_annot[[gene_ID_column]]),]


  pca_biplot = factoextra::fviz_pca_biplot(res.pca, geom = "text",
                               pointsize = 1 ,
                               labelsize = 2,
                               palette = colorsG,
                               col.var = gene_annot[[grouping_var]],
                               geom.ind = "point", col.ind = "lightgrey",
                               repel = TRUE,
                               select.var = list(contrib = 25),
                               show.legend = FALSE
  )+
    theme_classic(base_size = base_size)+
    labs(color=gene_ID_column)+
    theme_general()+
    theme(legend.position = "bottom",
          legend.direction = "vertical")+
    guides(color= guide_legend(
      direction = "vertical",
      title.position = "top",
      nrow=1, byrow=TRUE,
      title.theme = element_text(size = base_size, face = "bold"),
      label.theme = element_text(size = base_size)
    ))
  if (output == "biplot") {
    return(pca_biplot)
  }
  legend3 = cowplot::get_legend(pca_biplot)
  pca_both = cowplot::plot_grid(pca_ind, pca_biplot + theme(legend.position="none"))
  pca_full = cowplot::plot_grid(pca_both,
                                legend3,
                                rel_heights = c(6,1),
                                nrow = 2)
  return(pca_full)
}

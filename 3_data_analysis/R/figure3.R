



#' @title figure3
#' @description combines plots for figure 3
#'
#' @import cowplot
#' @export


figure3= function(){
  emptyggplot <- function() {
    return(ggplot() + theme_void())
  }
  PCAplt = plot_pca(expr_mat = exprmat_TMM_minmax,
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
                    output = "ind")+
    theme(legend.position = "none")
  clusterstats =  stats_E_cluster(embryo_annot = embryo_info,
                                  embryo_ID_column = "E_ID",
                                  selection_var = "Cells",
                                  color_var = "Morphology",
                                  fill_var= "ClusterNo",
                                  clusterEmbryos = clustering$Embryos,
                                  ClusterIDcolumnE = "ClusterNo",
                                  log_axis = T,
                                  colors = colorscheme,
                                  base_size = 8,
                                  horizontal_legend = F,
                                  xlab_txt = "Cluster No.",
                                  ylab_txt =  "Cells [#]",
                                  #"Days from Aug 1"^"st"
                                  color_lab_txt = "Morphology",
                                  fill_lab_txt = "Embryo\nCluster No."  ,
                                  grouponY = T)
  legendtop = get_legend(clusterstats+
                     theme(legend.position = "left",
                           legend.box="horizontal",
                           legend.box.just = "left"))
  clusterstats = clusterstats+theme(legend.position = "none")

  top_plots = plot_grid(plotlist = list(PCAplt,
                            emptyggplot(),
                            legendtop,
                            clusterstats),
            rel_widths = c(6.5,2,3.5,5.5),
            nrow = 1,
            labels = c("A", "", "", "", "B"),
            label_size = 12)

  expr_plots = expr_pattern_selected(expr_mat = exprmat_TMM_minmax,
                        gene_as_rows = T,
                        embryo_annot = embryo_info,
                        embryo_ID_column = "E_ID",
                        embryo_plot_var = "Cells",
                        gene_annot = markergenes,
                        gene_ID_column = "GeneName",
                        grouping_var = "show_pattern",
                        sorting_var = "plot_no",
                        color_var = "Compartment",
                        base_size = 8,
                        colors = colorscheme,
                        ylab_txt = "scaled TMM",
                        xlab_txt = "Cells [#]")

  expr_plots[["morph"]] = epiblast_morphology(ki_data = ki67evaluation,
                                              filter_col = "Stage",
                                              filter_var = "Diapausing",
                                              grouping_var = "ICM_shape",
                                              embryo_annot = embryo_info,
                                              x_var ="Nuclei_Tot",
                                              y_var = "Days_since_1.8.",
                                              colors = colorscheme[colorscheme$Class == "ICM",],
                                              base_size =8)
  center_plots = plot_grid(plotlist = expr_plots,
                                  labels = LETTERS[3:8], label_size = 12)
  micro_pics = assign(paste0("subfig"), image_read_svg(file.path(file.path("extdata", "images"), paste0("Fig3J.","svg")))) |> image_ggplot()

  combined = plot_grid(plotlist = list(top_plots, center_plots, micro_pics),
                       nrow = 3,
                       rel_heights = c(5.5,11,4.5))
  combined
  return(combined)

}





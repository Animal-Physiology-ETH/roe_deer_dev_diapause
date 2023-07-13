


#' @title stats_E_cluster
#' @description violin plots of embryo clusters (cells, days)
#'
#' @param embryo_annot data frame containing annotations for embryo data
#' @param embryo_ID_column column name with embryo / sample ID variable (must be identical to either column or row names of expr_mat)
#' @param selection_var variable in embryo_annot to analyze (numeric only)
#' @param color_var variable in embryo_annot to use for coloring points
#' @param fill_var variable in embryo_annot or clusterEmbryos to use for filling violoin plots
#'
#' @param clusterEmbryos data frame with cluster assignment for embryos / samples (ignored if compute_clusters_new = T)
#' @param ClusterIDcolumnE column name in clusterEmbryos with names for embryo / sample clusters
#'
#' @param log_axis plot on logarithmized axis
#' @param colors dataframe with color assignement
#' @param base_size base size for text (pt)
#' @param horizontal_legend plot legend horizontally
#' @param ylab_txt text on y axis
#' @param xlab_txt text on x axis
#' @param color_lab_txt text for color legend
#' @param fill_lab_txt text for fill legend
#' @param grouponY plot group on y axis
#'
#' @import dplyr
#' @export

stats_E_cluster = function(
    embryo_annot = embryo_info,
    embryo_ID_column = "E_ID",
    selection_var = "Cells",
    color_var = "Morphology",
    fill_var= "ClusterNo",

    clusterEmbryos = clustering$Embryos,
    ClusterIDcolumnE = "ClusterNo",

    log_axis = T,
    colors = colorscheme,
    base_size = 8,

    horizontal_legend = TRUE,
    xlab_txt = "Embryo Cluster No.",
    ylab_txt =  "Cells [#]",
    #"Days from Aug 1"^"st"
    color_lab_txt = "Morphology",
    fill_lab_txt = "Embryo\nCluster No."  ,
    grouponY = T){

  all_info = merge(embryo_annot, clusterEmbryos, by = embryo_ID_column)

  selection = unlist(lapply(all_info , is.numeric))
  selection[[embryo_ID_column]]=T
  selection[[ClusterIDcolumnE]]=T
  selection[[color_var]]=T

  all_numeric = all_info[,selection]

  summary_df = all_numeric %>%
    group_by_at(vars(one_of(ClusterIDcolumnE))) %>%
    summarize_at(vars(one_of(selection_var)),  list(mean = mean, sd = sd))

  if (horizontal_legend == TRUE) {
    guides = guides(fill = guide_legend(override.aes = list(shape = 16, size = 1),
                                        direction = "horizontal", title.position = "top",
                                        nrow=1, byrow=TRUE, order = 1,
                                        title.theme = element_text(size = base_size, face = "bold")),
                    color = guide_legend(override.aes = list(shape = 16, size = 1),
                                         direction = "horizontal", title.position = "top",
                                         nrow=1, byrow=TRUE, order = 2,
                                         title.theme = element_text(size = base_size, face = "bold")))
  }else{
    guides = guides(fill = guide_legend(override.aes = list(shape = 16, size = 1),
                                        direction = "vertical", title.position = "top",
                                        ncol = 1, byrow=T, order = 1,
                                        title.theme = element_text(size = base_size, face = "bold")),
                    color = guide_legend(override.aes = list(shape = 16, size = 1),
                                         direction = "vertical", title.position = "top",
                                         ncol = 1, byrow=T,order = 2,
                                         title.theme = element_text(size = base_size, face = "bold")))
  }
    plt = ggplot(data=all_numeric, aes_string(x=ClusterIDcolumnE,
                                                                  y=selection_var, group=ClusterIDcolumnE))+
      geom_violin(aes_string(fill = ClusterIDcolumnE), alpha = 0.8, width = 1)+
      geom_boxplot(width = 0.1, fill= "lightgrey")+
      geom_jitter(aes_string(colour = color_var), position=position_jitter(0.3), size = 1)+
      scale_color_manual(values = colors[colors[["Class"]]==color_var,][["Color"]])+
      scale_fill_manual(values = colors[colors[["Class"]]=="EmbryoCluster",][["Color"]])+
      xlab(xlab_txt) +
      theme_general()+
      theme(legend.position = "bottom",
            legend.direction = "vertical",
            legend.background=element_blank())+
      labs(color = color_lab_txt, fill = fill_lab_txt)+
      guides+
      ylab(ylab_txt)+
      add_pvalue(p_value_df(df = all_numeric,
                            variable = selection_var,
                            grouping = ClusterIDcolumnE,
                            log_scale = TRUE),
                 xmin = "group1",
                 xmax = "group2",
                 label = "stars",
                 y.position = "y.pos_remns",
                 cord.flip = grouponY)

    if (grouponY == T) {
        plt = plt+ coord_flip()
    }
    if (log_axis == T) {
      plt = plt + scale_y_continuous(trans = "log10",
                         breaks = scales::trans_breaks("log10",
                                                       function(x) 10^x),
                         labels = scales::trans_format("log10",
                                                       scales::math_format(10^.x)))+

      annotation_logticks(sides = ifelse(grouponY, "b", "l"))
    }
    return(plt)
}

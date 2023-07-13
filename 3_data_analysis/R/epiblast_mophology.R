#' @title epiblast_morphology
#' @description plots epiblast morphology and stat-boxes
#' @param ki_data data frame to use
#' @param filter_col column in ki_data to use for filtering data
#' @param filter_var variable(s) to filter ki_data for
#' @param grouping_var column in ki_data containing the groups
#' @param embryo_annot data frame containing annotations for embryo data
#' @param x_var variable in ki_data to plot on x axis
#' @param y_var variable in ki_data to plot on y axis
#' @param colors dataframe with color assignement
#' @param base_size base size for text (pt)
#' @returns ggplot object
#' @import ggplot2
#' @import dplyr
#' @export

epiblast_morphology = function(ki_data = ki67evaluation,
                               filter_col = "Stage",
                               filter_var = "Diapausing",
                               grouping_var = "ICM_shape",
                               embryo_annot = embryo_info,
                               x_var ="Nuclei_Tot",
                               y_var = "Days_since_1.8.",
                               colors = colorscheme[colorscheme$Class == "ICM",],
                               base_size =8
){
  #filter out elongated
  ki_data = ki_data[ki_data[[filter_col]] %in% filter_var,]
  groups = levels(factor(ki_data[[grouping_var]]))
  groups = groups[c(4,3,1,2)]
  ki_data[[grouping_var]] = factor(ki_data[[grouping_var]], levels = groups)

  summary_shape = ki_data %>%
    group_by_at(vars(one_of(grouping_var))) %>%
    summarize_at(vars(one_of(x_var, y_var)),  list(mean = mean, sd = sd))
  rownames(summary_shape) = summary_shape[[grouping_var]]

  min_days = summary_shape[which(summary_shape[[grouping_var]] =="round"),][[paste0(y_var, "_","mean")]]-
    summary_shape[which(summary_shape[[grouping_var]] =="round"),][[paste0(y_var, "_","sd")]]
  max_days = summary_shape[which(summary_shape[[grouping_var]] =="disk"),][[paste0(y_var, "_","mean")]]+
    summary_shape[which(summary_shape[[grouping_var]] =="disk"),][[paste0(y_var, "_","sd")]]

  p_dens = ggplot(data = ki_data,
                  aes_string(y = y_var, x=x_var ))+
    geom_point(aes_string(colour=grouping_var), size =2)+
    ylab("Days from Aug 1st")+
    xlab("Cells")+
    scale_color_manual(values = colors[["Color"]])+
    annotation_logticks(sides = "b")+
    scale_x_continuous(limits = c(90,
                                  max(embryo_annot$Cells) +
                                    0.1* max(embryo_annot$Cells)),
                       trans = "log10",
                       breaks = scales::trans_breaks("log10", function(x) 10^x),
                       labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_y_continuous(limits = c(min_days, max_days))+
    labs(colour="Epiblast\nMorphology")+
    geom_vline(xintercept =1000, linetype="dashed",
               color = "orange")+
    geom_vline(xintercept =7000, linetype="dashed",
               color = "orange")+
    theme_general()
  counter=0
  for (i in groups) {
    counter=+1
    sub_df = ki_data[ki_data[[grouping_var]] ==i,]
    n = length(sub_df)
    #for cells
    s_nuc = sd(sub_df[[x_var]])
    a_nuc = mean(sub_df[[x_var]])
    left_nuc = a_nuc-s_nuc
    rigth_nuc = a_nuc+s_nuc
    #for days
    s_d = sd(sub_df[[y_var]])
    a_d = mean(sub_df[[y_var]])
    left_d = a_d-s_d
    rigth_d = a_d+s_d

    p_dens = p_dens +
      annotate("rect", xmin = left_nuc,
               xmax = rigth_nuc,
               ymin = left_d, ymax = rigth_d,
               alpha = .2,fill = colors[colors[["Item"]]==i,][["Color"]])
  }
  p_dens = p_dens + guides( color = guide_legend(override.aes = list(shape = 16, size = 1),
                                                 direction = "vertical", title.position = "top",
                                                 ncol=1, byrow=TRUE,
                                                 label.theme = element_text(size = base_size),
                                                 title.theme = element_text(size = base_size, face = "bold")))
  p_dens = p_dens + theme(legend.position= c(0.6,0.8),
                          legend.spacing.y = unit(0.01, 'cm'))
  return(p_dens)
}

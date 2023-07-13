#' @title ki67_barplot
#' @description
#' plots barplot of ki67 counting data
#' @param ki_data data frame to use
#' @param grouping column of ki_data containing grouping variable
#' @param value_column column of ki_data containing value to display
#' @param ylabel label of y axis
#' @param colors data frame with colorscheme
#' @param colorclass class in colors to use
#' @importFrom scales trans_breaks
#' @importFrom scales trans_format
#' @importFrom scales math_format
#' @import ggplot2
#' @import svglite
#' @import dplyr
#' @importFrom ggprism add_pvalue
#' @export

ki67_barplot = function(ki_data = ki67evaluation,
                        grouping = "Stage",
                        value_column = "perc_tot_ki",
                        ylabel = "% of Ki67 positive Nuclei",
                        colors = colorscheme,
                        colorclass = "Morphology"){
  combined_df = as.data.frame(ki_data)
  combined_df$groups = combined_df[[grouping]]
  combined_df$value = combined_df[[value_column]]

  summary_df = combined_df %>%
    group_by(groups) %>%
    summarise(Ki67.positive = mean(100*value),
              Stdev = sd(100*value),
              SEM =std_error(100*value),
              n = n() )


  result = t.test(perc_tot_ki ~ groups, data = combined_df, var.equal = FALSE)$p.value
  result = signif(result, digits = 3)

  df_p_val = data.frame(
    group1 = unique(combined_df$groups)[1],
    group2 = unique(combined_df$groups)[2],
    label = result,
    label.stars = star_p(result),
    y.position = 85
  )
  coloringv = (colorscheme %>% filter(Class== colorclass))$Color

  bar_plotKi67 = ggplot(summary_df, aes(y=Ki67.positive, x=groups))+
    geom_errorbar(aes(ymin = Ki67.positive-Stdev,
                      ymax = Ki67.positive+Stdev),
                  width = 0.5,
                  size=0.5)+
    geom_col(aes(fill= groups), width = 0.8)+
    scale_y_continuous(expand = c(0,0)) +
    geom_text(aes(x=1, y=100.3, label="Stretch it"), vjust=-1)+
    geom_text(aes(y= min(Ki67.positive)/2,
                  label= paste0("n = ", n)),
              size = 3, col = "white",
              fontface = "bold")+
    scale_fill_manual(values = coloringv)+
    ylab(ylabel) +
    xlab(expression(paste("")))+
    labs(fill = "")+
    theme_general()+
    theme(legend.position = "none",
          legend.background=element_blank())

  bar_plotKi67 = bar_plotKi67+ ggprism::add_pvalue(df_p_val,
                                           xmin = "group1",
                                           xmax = "group2",
                                           label = "label.stars",
                                           y.position = "y.position",
                                           bracket.size = 0.4,
                                           label.size = 3.5,
  )
  return(bar_plotKi67)
}


#' @title ki67_box
#' @description
#' plots boxplot of ki67 counting data
#' @param ki_data data frame to use
#' @param grouping column of ki_data containing grouping variable
#' @param value_column column of ki_data containing value to display
#' @param ylabel label of y axis
#' @param colors data frame with colorscheme
#' @param colorclass class in colors to use
#' @importFrom scales trans_breaks
#' @importFrom scales trans_format
#' @importFrom scales math_format
#' @import ggplot2
#' @import svglite
#' @import dplyr
#' @importFrom ggprism add_pvalue
#' @export

ki67_box = function(ki_data = ki67evaluation,
                        grouping = "Stage",
                        value_column = "perc_tot_ki",
                        ylabel = "% of Ki67 positive Nuclei",
                        colors = colorscheme,
                        colorclass = "Morphology"){
  combined_df = as.data.frame(ki_data)
  combined_df$groups = combined_df[[grouping]]
  combined_df$value = combined_df[[value_column]]*100
  coloringv = (colorscheme %>% filter(Class== colorclass))$Color
  df_p_val = data.frame(
    group1 = unique(combined_df$groups)[1],
    group2 = unique(combined_df$groups)[2],
    label = result,
    label.stars = star_p(result),
    y.position = 75
  )
  ki67_box = ggplot(combined_df, aes(x = groups, y = value))+
    geom_boxplot(aes(fill = groups))+
    scale_y_continuous(expand = c(0,1.2)) +
    scale_fill_manual(values = coloringv)+
    ylab(ylabel) +
    xlab(expression(paste("")))+
    labs(fill = "")+
    theme_general()+
    theme(legend.position = "none",
          legend.background=element_blank())


  ki67_box = ki67_box+ ggprism::add_pvalue(df_p_val,
                                                   xmin = "group1",
                                                   xmax = "group2",
                                                   label = "label.stars",
                                                   y.position = "y.position",
                                                   bracket.size = 0.6,
                                                   label.size = 4,
  )
}




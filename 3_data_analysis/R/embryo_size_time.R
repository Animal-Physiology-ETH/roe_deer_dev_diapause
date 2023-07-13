#' @title embryo_cell_time
#' @description plots the embryo cells against time
#' @param e_data data frame to use
#' @param xparam column to plot on x axis
#' @param coloring column to use for coloring
#' @param cells_on_x display cells on x axis (T/F) if == T plot is flipped
#' @param colors data frame with colorscheme
#' @param colorclass class in colors to use
#' @param xlabel label of x axis
#' @param legend_pos position of legend (vector of lenght 2, with values from 0 to 1)
#' @importFrom scales trans_breaks
#' @importFrom scales trans_format
#' @importFrom scales math_format
#' @import ggplot2
#' @export

embryo_cell_time = function(e_data = embryo_info,
                            xparam = "Days_since_1.8.",
                            coloring = "Morphology",
                            cells_on_x = F,
                            colors = colorscheme,
                            colorclass = "Morphology",
                            xlabel = bquote('Days from '~1^"st"~"of August"),
                            legend_pos = c(0.2, 0.9)
                            ){
  coloringv = (colorscheme %>% filter(Class== colorclass))$Color
  e_data$x = e_data[[xparam]]
  e_data$color =e_data[[coloring]]

  cellno_vs_date <-
    ggplot(e_data,aes(y=Cells, x=x, colour = color)) +
    geom_point(size = 2)+
    scale_y_continuous(trans = "log10",
                       limits = c(100,1000000),
                       breaks = scales::trans_breaks("log10", function(x) 10^x),
                       labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_color_manual(values = coloringv)+
    xlab(xlabel)+
    ylab("Estimated Cell Number")+
    labs(color = "")+
    theme_general()+
    theme(legend.position = legend_pos)
  if (cells_on_x ==T) {
    cellno_vs_date = cellno_vs_date + coord_flip()+annotation_logticks(sides = "b")
  }else{
    cellno_vs_date = cellno_vs_date + annotation_logticks(sides = "l")
  }
  return(cellno_vs_date)
}


#' @title embryo_cell_time_alt
#' @description plots the embryo cells against time (alternative with including cluster)
#' @param e_data data frame to use
#' @param xparam column to plot on x axis
#' @param coloring column to use for coloring
#' @param cells_on_x display cells on x axis (T/F) if == T plot is flipped
#' @param colors data frame with colorscheme
#' @param colorclass class in colors to use
#' @param xlabel label of x axis
#' @param legend_pos position of legend (vector of lenght 2, with values from 0 to 1)
#' @importFrom scales trans_breaks
#' @importFrom scales trans_format
#' @importFrom scales math_format
#' @import ggplot2
#' @export

embryo_cell_time_alt = function(e_data = embryo_info,
                            xparam = "Days_since_1.8.",
                            coloring = "ClusterNo",
                            shapes = "Morphology",
                            cells_on_x = F,
                            colors = colorscheme,
                            colorclass = "EmbryoCluster",
                            xlabel = bquote('Days from '~1^"st"~"of August"),
                            legend_pos = c(0.05, 0.95)
){
  coloringv = (colorscheme %>% filter(Class== colorclass))$Color
  e_data = merge(e_data, clustering$Embryos, by.y = "E_ID", by.x = "E_ID", all.x = T)
  e_data$x = e_data[[xparam]]
  e_data$color =e_data[[coloring]]
  e_data$shape =e_data[[shapes]]

  cellno_vs_date <-
    ggplot(e_data,aes(y=Cells, x=x, colour = color, shape = shape)) +
    geom_point(size = 2)+
    scale_y_continuous(trans = "log10",
                       limits = c(100,1000000),
                       breaks = scales::trans_breaks("log10", function(x) 10^x),
                       labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_color_manual(values = coloringv)+
    xlab(xlabel)+
    ylab("Estimated Cell Number")+
    labs(shape = "", color = "Cluster #")+
    theme_general()+
    theme(legend.position = legend_pos,
          legend.direction = "horizontal",
          legend.box.background = element_blank())
  if (cells_on_x ==T) {
    cellno_vs_date = cellno_vs_date + coord_flip()+annotation_logticks(sides = "b")
  }else{
    cellno_vs_date = cellno_vs_date + annotation_logticks(sides = "l")
  }
  return(cellno_vs_date)
}

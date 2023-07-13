#' @title theme_general
#' @description
#' general theme for ggplots
#' @param textsize base size of text
#' @import ggplot2
#' @export

theme_general <- function(textsize = 8){
  theme_minimal() %+replace%
    theme(legend.background = element_blank(),
          legend.box.background = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.justification = c("left", "top"),
          legend.key=element_blank(),
          text=element_text(size=textsize),
          axis.text=element_text(size=textsize),
          axis.title=element_text(size=textsize),
          plot.title=element_text(size=textsize+1),
          legend.text=element_text(size=textsize),
          legend.title=element_text(size=textsize),
          legend.key.size = unit(0.7, 'lines'))
}

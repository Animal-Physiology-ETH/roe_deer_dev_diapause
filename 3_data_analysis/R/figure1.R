##### 1 #####
#' @title figure1
#' @description generates figure 1 for the manuscript
#' @param image_dir directory with images
#' @param h plot height
#' @param w plot width
#' @param fontsize_labels size of labels (pt)
#' @param save save plot? (T/F)
#' @returns ggplot object
#' @import ggplot2
#' @import cowplot
#' @import magick
#' @export

figure1 = function(image_dir = file.path("extdata", "images"),
                       fontsize_labels = 12) {
  for (subfig in c("B","C","D","F")) {
    assign(paste0("F1", subfig), image_read(file.path(image_dir, paste0("Fig1", subfig,".png")), density = 600) |> image_ggplot())
  }

  F1F_ = cowplot::plot_grid(F1F, labels = "F", label_colour = "white", label_size = fontsize_labels, label_x = 0.01)

  plots_right =  cowplot::plot_grid(F1D,
                                    ki67_box(),
                                    F1F_,
                                    ncol = 1,
                                    rel_heights = c(45,60,45),
                                    labels = c("D", "E"),
                                    label_size = fontsize_labels,
                                    label_x = 0.01)

  plots_center = cowplot::plot_grid(F1B, ggplot()+theme_void(),
                                    F1C, ggplot()+theme_void(),
                                    ncol = 1,
                                    labels  = c("B", "", "C", ""),
                                    rel_heights = c(95,15,35,10),
                                    label_size = fontsize_labels)

  F1_complete = cowplot::plot_grid(embryo_cell_time(),plots_center,  plots_right, nrow = 1,
                                   rel_widths = c(100,35, 45),
                                   labels  = c("A", "", ""),
                                   label_size = fontsize_labels)
  return(F1_complete)
}

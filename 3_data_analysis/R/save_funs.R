#' @title save_heatmap
#' @description
#' draws heatmap with the option to save it.
#' @param h height in cm
#' @param w width in cm
#' @param outdir output directory
#' @param outformat format for output (svg, pdf, png)
#' @param filename add custom filename
#' @param ht heatmap to print or save
#' @param resolution resolution (dpi)
#' @param seed seed used for kmeans clustering
#' @importFrom svglite svglite
#' @import grid
#' @import ComplexHeatmap
#' @export

save_heatmap = function(w = 18,
                        h =  22,
                        outdir = "output",
                        filename = "Heatmap",
                        outformat = "svg",
                        ht = custom_heatmap(),
                        resolution = 600,
                        seed = 123){
  file = file.path(outdir, paste0(filename, ".", outformat))
  if (file.exists(file)) {
    file = file.path(outdir,
                     paste0(paste0(filename,"_", format(Sys.Date(), '%y%m%d'),
                                   "_", format(Sys.time(), "%H%M")),
                            ".", outformat))

  }
  if (outformat == "svg") {
    svglite::svglite(filename = file, width = w/2.54, height =  h/2.54)
    set.seed(seed)
    draw(ht,  merge_legends = TRUE)
    dev.off()
  }
  if (outformat == "pdf") {
    cairo_pdf(filename = file, width = w/2.54, height =  h/2.54, fallback_resolution = resolution)
    set.seed(seed)
    draw(ht,  merge_legends = TRUE)
    dev.off()
  }
  if (outformat == "png") {
    png(filename = file, width = cmtopixel(w), height =  cmtopixel(h), res = resolution)
    set.seed(seed)
    draw(ht,  merge_legends = TRUE)
    dev.off()
  }
  if (outformat == "tif") {
    tiff(filename = file, width = cmtopixel(w, resolution),
         height =  cmtopixel(h, resolution), res = resolution,
         compression = "none")
    set.seed(seed)
    draw(ht,  merge_legends = TRUE)
    dev.off()
    } else{
    set.seed(seed)
    draw(ht,  merge_legends = TRUE)
  }
}


#' @title save_plt
#' @description
#' plots and saves any plot
#' @param plt plot to save
#' @param w width in cm
#' @param h height in cm
#' @param outdir output directory
#' @param outformat format for output (svg, pdf, png)
#' @param filename add custom filename
#' @param resolution resolution (dpi)
#' @importFrom svglite svglite
#' @import grid
#' @import ggplot2
#' @export

save_plt = function(plt = ggplot(),
                    w = 18,
                    h =  22,
                    outdir = "output",
                    filename = "",
                    outformat = "pdf",
                    resolution = 600){
  file = file.path(outdir, paste0(filename, ".", outformat))
  if (file.exists(file)) {
    file = file.path(outdir,
                     paste0(paste0(filename,"_", format(Sys.Date(), '%y%m%d'),
                                   "_", format(Sys.time(), "%H%M")),
                            ".", outformat))

  }
  if (outformat == "svg") {
    svglite::svglite(filename = file, width = w/2.54, height =  h/2.54)
    plot(plt)
    dev.off()

  }
  if (outformat == "pdf") {
    cairo_pdf(filename = file, width = w/2.54, height =  h/2.54, fallback_resolution = resolution)
    plot(plt)
    dev.off()
  }
  if (outformat == "png") {
    png(filename = file, width = cmtopixel(w, resolution), height =  cmtopixel(h, resolution), res = resolution)
    plot(plt)
    dev.off()
  }
  if (outformat == "tif") {
    tiff(filename = file, width = cmtopixel(w, resolution),
         height =  cmtopixel(h, resolution), res = resolution,
         compression = "none")
    plot(plt)
    dev.off()
  } else{
    plot(plt)
  }
}


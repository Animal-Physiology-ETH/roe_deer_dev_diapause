
#' @title cmtopixel
#' @description
#' converts cm to pixelvalues at res dpi
#' @param cm cm to convert
#' @param res resolution
#' @export

cmtopixel = function(cm, res = 600){return((res * cm)/2.54)}


#' @title cmtoinch
#' @description
#' converts cm to inch
#' @param cm cm to convert
#' @export

cmtoinch = function(cm){return((cm)/2.54)}

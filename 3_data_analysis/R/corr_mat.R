#' @title mycolor_label
#' @description
#' creates custom color label legend
#' @param color_annotation color annotation df
#' @param label_size relative size of labels
#' @param point_size length of the sides of the squares
#' @param x_shift shift of legend in x direction
#' @param y_shift shift of legend in y direction
#' @param y_distance distance between legend lines

mycolor_label <- function (color_annotation,
                           label_size = .8,
                           point_size = 0.5,
                           x_shift = -6,
                           y_shift = 2,
                           y_distance = 1.5) {
  for (lbl in names(color_annotation)) {
    symbols(x_shift, y_shift, squares = point_size,
            add = TRUE, inches = FALSE,
            bg = color_annotation[lbl],
            fg = color_annotation[lbl])
    text((x_shift+2*label_size), (y_shift+0.5), labels = lbl,
         adj = c(0,1), cex = label_size)
    y_shift = y_shift +y_distance
  }
}

#' @title light_dark
#' @description
#' returns black for light colors and white for dark colrs to allow better readability
#' @param hexcolor color in hex

light_dark = function(hexcolor){
  #convert to RGB
  r = col2rgb(hexcolor)[1]
  g = col2rgb(hexcolor)[2]
  b = col2rgb(hexcolor)[3]
  #HSP equation from http://alienryderflex.com/hsp.html
  hsp = sqrt(
    0.299 * (r * r) +
      0.587 * (g * g) +
      0.114 * (b * b)
  )
  #Using the HSP value, determine whether the color is light or dark
  if (hsp>127.5) {
    return("#000000")
  }else {
    return("#FFFFFF")
  }
}


#' @title blocks_names
#' @description
#' creates boxes around text labels with colors indicated in named_colorvector
#' @param named_colorvector named vector containing colors
#' @param border_color border color of boxes
#' @param label_size relative size of labels

blocks_names <- function(named_colorvector,
                         border_color = NA,
                         label_size = .4) {
  drawcell <- function(fx, fy, bgcolor, text1, txtcol) {

    oldpar <- par(mar = c(0,0,0,0), bg = "white")
    on.exit(par(oldpar), add = TRUE)
    text(fx, fy, labels = text1, col = txtcol, cex = label_size)
  }

  drawcell_big <- function(fx, fy, bgcolor, text1, txtcol) {

    oldpar <- par(mar = c(0,0,0,0), bg = NA)
    on.exit(par(oldpar), add = TRUE)
    symbols(fx, fy, rectangles = matrix(c(7,1), 1, 2), #matrix(c(7,1), nrow = 1, ncol = 2),
            add = TRUE, inches = FALSE, fg = border_color, bg = bgcolor)
  }

  for (x in 1:length(named_colorvector)) {
    for (y in x:length(named_colorvector)) {
      if (x == y) {
        drawcell_big(x-3, length(named_colorvector) - y + 1, bgcolor = named_colorvector[x], text1 =NA, txtcol = NA)
        drawcell(x-3, length(named_colorvector) - y + 1, bgcolor = NA, text1 = names(named_colorvector)[x], txtcol = light_dark(named_colorvector[x]))
      }
    }
  }
}





#' @title corrmat
#' @description
#' Plots expression levels of selected df
#' @param expr_df vector of selected genes
#' @param gene_annot data frame containing annotations for genes (also defines which genes to display)
#' @param gene_ID_column column name with gene ID (must be identical to either column or row names of expr_mat).
#' @param colors dataframe with color assignement
#' @param colors_chosen refers to Item column in colors to be display (Process or Compartment)
#' @param label_size size of labels relative to default
#' @import corrplot
#' @import RColorBrewer
#' @importFrom Hmisc rcorr
#' @export

corrmat <- function(expr_mat = exprmat_TMM_minmax,
                    gene_annot = markergenes,
                    gene_ID_column = "GeneName",
                    colors = colorscheme,
                    colors_chosen = "Compartment",
                    label_size = 0.4){
  gene_annot = gene_annot[gene_annot[[gene_ID_column]] %in% rownames(expr_mat),]
  expr_mat = expr_mat[rownames(expr_mat) %in% gene_annot[[gene_ID_column]],]
  expr_mat <- t(expr_mat[order(rownames(expr_mat)) , ])
  order <- corrMatOrder(rcorr(expr_mat, type = "spearman")$r, order = "hclust")

  #load annotation df
  color_df = merge(gene_annot, colors, by.x = colors_chosen, by.y = "Item")
  color_df = color_df[order(color_df[[gene_ID_column]]),]
  color_df = color_df[order,]


  color_vector_genes = color_df[["Color"]]
  color_vector_groups = color_vector_genes
  names(color_vector_genes) = color_df[[gene_ID_column]]

  color_vector_groups = unique(color_df[order(color_df[["Compartment_No"]], decreasing = T),c("Color",colors_chosen)])
  color_vector_groups = setNames(color_vector_groups[["Color"]], color_vector_groups[[colors_chosen]])

  rhos <- rcorr(as.matrix(expr_mat), type = "spearman")$r
  rhos[is.na(rhos)] <- 0
  pval <- rcorr(as.matrix(expr_mat), type = "spearman")$P
  pval[is.na(pval)] <- 1
  corrplot(rhos,
             mar = c(5, 5,5, 5),
             p.mat = pval,
             order = "hclust",
             hclust.method = "complete",
             #bg = "blue",
             tl.col = color_vector_genes,
             insig = "blank", pch.col = "black",
             tl.srt = 90,
             pch.cex = label_size,
             sig.level = c(.05),
             tl.cex = label_size,
             cl.cex = label_size,
             cl.ratio = 0.1,
           shade.col	=NA,

             col  = colorRampPalette(brewer.pal(11, "Spectral"))(100),
             type = "upper")
    mycolor_label(color_vector_groups,  x_shift = 5, y_shift = 15, label_size = label_size)
    blocks_names(color_vector_genes, label_size = label_size)
}

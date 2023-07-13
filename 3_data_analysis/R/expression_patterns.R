
#' @title expr_vs_cells
#' @description
#' plots gene expression vs embryo cell number
#' @param expr_mat numeric matrix containing expression data, colnames (or rownames) should contain values found in embryo_ID_column
#' @param gene_as_rows are genes as rows? (T/F) --> if false matrix will by transposed
#'
#' @param embryo_annot data frame containing annotations for embryo data
#' @param embryo_plot_var column in embryo annot to plot on x
#'
#' @param gene_annot data frame containing annotations for genes (also defines which genes to display)
#' @param  vertical_lines display vertical lines as orientation
#' @param points display individual data points

#' @param ylab_txt text on y axis
#' @param xlab_txt text on x axis

#' @param linecol color for lines
#' @param filling = color for shades (CI)

#' @param base_size base size for text (pt)

#' @import ggplot2
#' @importFrom scales math_format

expr_vs_cells <- function(gene = "NANOG",
                          expr_mat = exprmat_TMM_minmax,

                          embryo_annot = embryo_info,
                          embryo_plot_var = "Cells",

                          gene_annot = markergenes,

                          base_size = 8,

                          vertical_lines = TRUE,
                          points = FALSE,

                          ylab_txt = "scaled TMM",
                          xlab_txt = "Cells [#]",

                          linecol = "black",
                          filling = "black"
                          ){
  expr_df = as.data.frame(t(expr_mat))
  p <- ggplot(expr_df, aes_string(x=embryo_annot[[embryo_plot_var]], y = gene))+
    theme_general()+
    geom_smooth(stat="smooth",method = "loess", formula= y ~ x, color = linecol, alpha = 0.2, fill = filling)+
    scale_x_continuous(limits = c(min(embryo_annot[[embryo_plot_var]])/2,
                                  1.1* max(embryo_annot[[embryo_plot_var]])),
                       trans = "log10",
                       breaks = scales::trans_breaks("log10", function(x) 10^x),
                       labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    ylab(ylab_txt)+
    xlab(xlab_txt)+
    labs(color = "")+
    annotation_logticks(sides = "b")

  if (points == TRUE) {
    p = p +geom_point()
  }

  if (vertical_lines == TRUE) {
    p = p +
      geom_vline(xintercept =1000, linetype="dashed",
                 color = "orange")+
      geom_vline(xintercept =7000, linetype="dashed",
                 color = "orange")
  }

  return(p)
}


#' @title expr_pattern_selected
#' @param expr_mat numeric matrix containing expression data, colnames (or rownames) should contain values found in embryo_ID_column
#' @param gene_as_rows are genes as rows? (T/F) --> if false matrix will by transposed
#'
#' @param embryo_annot data frame containing annotations for embryo data
#' @param embryo_ID_column column name with embryo / sample ID variable (must be identical to either column or row names of expr_mat)
#' @param embryo_plot_var variable in embryo_annot to plot (numeric)
#'
#' @param gene_annot data frame containing annotations for genes (also defines which genes to display)
#' @param gene_ID_column column name with gene ID (must be identical to either column or row names of expr_mat).
#' Order of elements in gene_ID_column will be maintained inside each gene cluster.
#' @param grouping_var variable in gene_annot containing the groups to be displayed
#' @param sorting_var variable in gene_annot containing the order of plots to be displayed
#' @param color_var variable in gene_annot to use for coloring (must be "Class" in colors)

#' @param base_size base size for text (pt)
#' @param colors dataframe with color assignement
#' @param ylab_txt text on y axis
#' @param xlab_txt text on x axis
#'
#' @import factoextra
#' @import ggplot2
#'
expr_pattern_selected = function(expr_mat = exprmat_TMM_minmax,
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
                               xlab_txt = "Cells [#]"
                               ){
  if(gene_as_rows ==F){
    expr_mat=t(expr_mat)
  }
  colordf = merge(gene_annot, colors[colors$Class == color_var,], by.x = color_var, by.y = "Item")
  groups = unique(gene_annot[order(gene_annot[[sorting_var]]),][[grouping_var]])
  groups = groups[!is.na(groups)]

  plot_collection = list()

  for (group in groups) {
    d=0
    genes = gene_annot[which(gene_annot[[grouping_var]]==group),gene_ID_column][[1]]

    compartments = unique(gene_annot[gene_annot[[gene_ID_column]] %in% genes, color_var])[[1]]
    coloring = unique(colordf[colordf[[color_var]] %in% compartments, "Color"])

    colormix = coloring[1]
    if (length(coloring) >1) {
      for (c in 1:(length(coloring)-1)) {
        colormix = colorspace::hex(colorspace::mixcolor(0.5, colorspace::hex2RGB(coloring[c]),
                                                        colorspace::hex2RGB(coloring[c+1])))
        coloring[c+1]=colormix
      }
    }

    pl = expr_vs_cells(gene = genes[1],
                       expr_mat = expr_mat,
                       embryo_annot = embryo_annot,
                       embryo_plot_var = "Cells",
                       gene_annot = gene_annot,
                       base_size = base_size,
                       ylab_txt = "scaled TMM",
                       xlab_txt = "Cells [#]",
                       linecol = NA,
                       filling = colordf[colordf[[gene_ID_column]] == genes[1], "Color"])+
      annotate("text", x = 0.95*max(embryo_annot[[embryo_plot_var]]),
               y = max(1)-d,
               label = paste(genes[1]), vjust = 1,
               hjust = 1, fontface ="italic", size = 8/.pt,
               color = "black") #gene_annot[genes[1], "mycolors"]) #change here and within loop to get colors in names back!

    for (g in genes[2:length(genes)]) {
      d=d+0.1
      pl = pl + geom_smooth(aes_string(y = t(expr_mat)[,g], x = embryo_annot[[embryo_plot_var]]),
                            stat="smooth",method = "loess", formula= y ~ x,
                            color = NA,
                            fill = colordf[colordf[[gene_ID_column]] == genes[1], "Color"], alpha = 0.2,
                            size =  0.5)+
        annotate("text", x = 0.95*max(embryo_annot[[embryo_plot_var]]),
                 y = max(1)-d,
                 label = paste(g), vjust = 1,
                 hjust = 1, fontface ="italic", size = 8/.pt,
                 color = "black") #gene_annot[genes[1], "mycolors"]) #change here and above to get colors in names back!
    }
    expr_group <- expr_mat[rownames(expr_mat) %in% genes,]
    expr_group = colMeans(expr_group)
    pl = pl + geom_line(aes_string(y = expr_group, x = embryo_annot[[embryo_plot_var]]),
                        stat="smooth",method = "loess", formula= y ~ x, color = colormix, size = 1)
    pl = pl + ylab(ylab_txt)+xlab(xlab_txt)
    pl= pl +
      scale_y_continuous(limits = c(-0.2, 1.1))+
      theme_general()
    plot_collection[[group]] = pl
  }
  return(plot_collection)
  }

#' @title expression_swarm_plt
#' @description
#' Plots expression levels of selected df
#' @param expr_mat matrix with expression data (not scaled!)
#' @param gene_annot data frame containing annotations for genes (also defines which genes to display)
#' @param gene_ID_column column name with gene ID (must be identical to either column or row names of expr_mat).
#' @param colors dataframe with color assignement
#' @param colors_chosen refers to Item column in colors to be display (Process or Compartment)
#' @param sorting_col column to sort plot by
#' @param logscale do you want a log-y axis ? (T/F)
#' @param allin1 should all plots be in one ? (T/F) (if F will generate four plots via quantiles)
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom scales trans_format
#' @importFrom scales math_format
#' @importFrom ggbeeswarm geom_quasirandom
#' @export

expression_swarm_plt <- function(expr_mat = exprmat_TMM,
                                 gene_annot = markergenes,
                                 gene_ID_column = "GeneName",
                                 colors = colorscheme,
                                 colors_chosen = "Compartment",
                                 sorting_col = "Compartment_No",
                                 logscale = T,
                                 allin1 = F){

  gene_annot = gene_annot[gene_annot[[gene_ID_column]] %in% rownames(expr_mat),]
  expr_mat = expr_mat[rownames(expr_mat) %in% gene_annot[[gene_ID_column]],]
  expr_df <- as.data.frame(t(expr_mat[order(rownames(expr_mat)) , ]))

  reshape_df = function(expr_df = expr_df){
    expr_df.long<-reshape2::melt(t(expr_df),id.vars = colnames(expr_df))
    expr_df.long$value = as.numeric(expr_df.long$value)
    expr_df.long = merge(expr_df.long, gene_annot, by.x = "Var1", by.y = gene_ID_column, all.x = T)
    expr_df.long = merge(expr_df.long, colors, by.x = colors_chosen, by.y = "Item")
    expr_df.long$Color[which(is.na(expr_df.long$Color))] <- "#000000"

    expr_df.long[[colors_chosen]][which(is.na(expr_df.long[[colors_chosen]]))] <- "None"
    expr_df.long <- expr_df.long[order(expr_df.long[[sorting_col]]),]

    expr_df.long$Var1 = factor(expr_df.long$Var1,levels = as.factor(unique(expr_df.long[order(expr_df.long[[sorting_col]]),][["Var1"]])))

    expr_df.long$MarkerGroup2 <- expr_df.long[[colors_chosen]]
    return(expr_df.long)
  }

  expr_df.long = reshape_df(expr_df)
  color_v = setNames(unique(expr_df.long$Color), unique(expr_df.long$MarkerGroup2))



  if (allin1 == F) {
    quantiles = function(input_df = expr_df){
      max_count =  input_df %>% summarise_if(is.numeric, max)
      mean_count =  input_df %>% summarise_if(is.numeric, mean)
      min_count = input_df %>% summarise_if(is.numeric, min)
      summary_df = rbind("min" = min_count, "max" = max_count,"mean" = mean_count)
      summary_df["range", ] =  summary_df["max", ]- summary_df["min", ]
      summary_df = as.data.frame(t(summary_df))
      q = quantile(summary_df$mean)
      gene_quantiles = list(
        "superlow" = rownames(summary_df %>% filter(mean < q[2])),
        "low" = rownames( summary_df %>% filter((mean >= q[2]) & (mean < q[3]))),
        "med" = rownames(summary_df %>% filter((mean >= q[3]) & (mean < q[4]))),
        "high" = rownames(summary_df %>% filter(mean >= q[4]))
      )
      return(gene_quantiles)
    }
    gene_quantiles = quantiles()

    plt_ls = list()
    for (quantile in names(gene_quantiles)) {
      expr_sel = expr_df.long %>% filter(Var1 %in% gene_quantiles[[quantile]])
      expr_sel = expr_sel %>% filter(value !=0)
      exprlvl_swarm <-  ggplot(expr_sel, aes(y=value,
                                             x= Var1,
                                             group = Var1,
                                             color =MarkerGroup2))+
        geom_quasirandom()+
        ylab("TMM counts")+
        theme_general()+
        theme(axis.line = element_line(colour = "black"),
              axis.text.x = element_text(angle = 90, vjust = 0.5,
                                         hjust=1),
              legend.key=element_blank(),
              legend.position = "bottom",
              legend.box="vertical", legend.margin=margin(),
              axis.title.x = element_blank())+
        scale_color_manual(values = color_v)

      if (logscale == T) {
        exprlvl_swarm = exprlvl_swarm +
          scale_y_continuous(trans = "log10",
                             labels = scales::trans_format("log10", scales::math_format(10^.x)))
      }
      plt_ls[[quantile]] = exprlvl_swarm + theme(legend.position = "none")
    }
    plts = plot_grid(plotlist = plt_ls)
  }

  plt_all = ggplot(expr_df.long, aes(y=value,
                                 x= Var1,
                                 group=Var1,
                                 color =MarkerGroup2))+
    geom_jitter()+
    ylab("TMM counts")+
    theme_general()+
    theme(axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust=1),
          legend.key=element_blank(),
          legend.position = "bottom",
          legend.box="vertical", legend.margin=margin(),
          axis.title.x = element_blank())+
      scale_color_manual(values = color_v)+
    labs(color = "")
  legend=get_legend(plt_all)

  if (allin1 == T) {
    if (logscale == T) {
      plt_all = plt_all +
        scale_y_continuous(trans = "log10",
                           labels = scales::trans_format("log10", scales::math_format(10^.x)))
    }
    return(plt_all)
  }else{
    plts_incl = plot_grid(plotlist = list(plts, legend), ncol = 1, rel_heights = c(8,1))
    return(plts_incl)
  }

}

#4f6376


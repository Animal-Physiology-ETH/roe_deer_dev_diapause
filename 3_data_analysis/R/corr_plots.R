
#' @title corr_groups
#' @description
#' form groups with high correlations
#' @param expr_df dataframe with gene expression (columnnames must be genes)
#' @param corr_type type of correlation to use in rcorr (spearman, pearson, ...)
#' @param R_limit R/Rho threshold
#' @param p_limit p value threshold
#' @param h_cut height to cut the hierarchical clustering tree
#' @param n_group number of groups (NULL if group number should not be pre-determined)
#' @importFrom Hmisc rcorr
#' @export

corr_groups = function(expr_df,
                       corr_type = "spearman",
                       R_limit = 0.7,
                       p_limit = 0.05,
                       h_cut = 0.7,
                       n_group = NULL){
  rho = as.data.frame(rcorr(as.matrix(expr_df), type = corr_type)$r)
  pval = as.data.frame(rcorr(as.matrix(expr_df), type = corr_type)$P)
  rho = as.data.frame(rcorr(as.matrix(expr_df), type = corr_type)$r)
  pval = as.data.frame(rcorr(as.matrix(expr_df), type = "spearman")$P)
  rho_t = rho
  rho_t[(pval > p_limit)] <- 1
  rho_t[(abs(rho_t)< R_limit)] <- NA
  rho_t[rho_t == 1] <- NA

  rho_filtered = rho %>% select(colnames(rho_t %>% select_if(~ !all(is.na(.)))))
  expr_filtered = expr_df %>% select(colnames(rho_t %>% select_if(~ !all(is.na(.)))))

  set.seed(123)
  distance_mat = dist(t(rho_filtered), method = 'maximum')
  hier_clus = hclust(distance_mat, method = "complete")
  #plot(hier_clus)
  fit <- cutree(hier_clus , h = h_cut, k =n_group)
  fit_df = data.frame("GeneName" = names(fit), "ClusterNo" = fit)
  fit_df = fit_df %>% filter(GeneName %in% colnames(rho_filtered))
  return(fit_df)

}


#' @title correval_plts
#' @description
#' plots gene groups with high correlation
#' @param expr_mat matrix with gene expression (genes as rows)
#' @param gene_annot data frame containing annotations for genes (also defines which genes to display)
#' @param gene_ID_column column name with gene ID (must be identical to either column or row names of expr_mat).
#' @param corr_type type of correlation to use in rcorr (spearman, pearson, ...)
#' @param coef_name name of coefficient in plots
#' @param R_limit R/Rho threshold
#' @param p_limit p value threshold
#' @param h_cut height to cut the hierarchical clustering tree
#' @param n_group number of groups (NULL if group number should not be pre-determined)
#' @param labels vector with labels (AUTO/auto for big/small automatic letters)
#' @importFrom Hmisc rcorr
#' @export

correval_plts= function(expr_mat = exprmat_TMM_minmax,
                         gene_annot = markergenes,
                         gene_ID_column = "GeneName",
                         corr_type = "spearman",
                         coef_name = "rho",
                         R_limit = 0.7,
                         p_limit = 0.05,
                         h_cut = 0.7,
                         n_group = NULL,
                         labels = "AUTO"){

  gene_annot = gene_annot[gene_annot[[gene_ID_column]] %in% rownames(expr_mat),]
  expr_mat = expr_mat[rownames(expr_mat) %in% gene_annot[[gene_ID_column]],]
  expr_df <- as.data.frame(t(expr_mat[order(rownames(expr_mat)) , ]))

  plot_coll = list()
  fit_df = corr_groups(expr_df = expr_df,
                       corr_type = corr_type,
                       R_limit = R_limit,
                       p_limit = p_limit,
                       h_cut = h_cut,
                       n_group = n_group)
  fit_df$ClusterNo= LETTERS[fit_df$ClusterNo]
  fit_expr_df = merge(fit_df, t(expr_df), by= 0) %>%
    select(-c(Row.names))

  group_means = fit_expr_df %>%
    group_by(ClusterNo) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

  group_means_t = as.data.frame(t(group_means %>% select(-c(ClusterNo)) ))
  colnames(group_means_t) = paste0("CG", group_means$ClusterNo)


  rho = as.data.frame(rcorr(as.matrix(group_means_t), type = corr_type)$r)
  pval = as.data.frame(rcorr(as.matrix(group_means_t), type = corr_type)$P)

  plot_coll = list()
  for(i in LETTERS[1:(ncol(group_means_t)-1)]){
    for (j in LETTERS[2:(ncol(group_means_t))]) {
      if (i!=j && i < j) {
        rho_c = abs(rho[paste0("CG",i), paste0("CG",j)])
        pval_c = (pval[paste0("CG",i), paste0("CG",j)])
        if (rho_c >= R_limit && pval_c < p_limit) {
          namex = ifelse(j == "A", "naive group",
                         ifelse(j == "B", "mixed group",
                                ifelse(j == "C", "core/epiblast group",
                                       ifelse(j == "D", "ex. endoderm group",
                                              ifelse(j == "F", "trophectoderm group",""
                                              )))))
          namey = ifelse(i == "A", "naive group",
                         ifelse(i == "B", "mixed group",
                                ifelse(i == "C", "core/epiblast group",
                                       ifelse(i == "D", "ex. endoderm group",
                                              ifelse(i == "F", "trophectoderm group",""
                                              )))))
          pl = ggplot(group_means_t, aes_string(x = paste0("CG",j), y = paste0("CG",i)))+
            stat_smooth(size = 0.5,
                        method = "lm",
                        formula = y~x, fill = "lightgrey", color = "#d49494")+
            geom_point(size = 0.5)+
            theme_general()+
            theme(plot.title = element_text(face = "bold"))+
            ggpubr::stat_cor(aes(label =  paste(..r.label..)),
                             method = corr_type, cor.coef.name = coef_name ,
                             label.x.npc = 1, label.y.npc = 1, hjust = 1, vjust = 1)+
            #labs(title = paste("Correlation group", i, "vs", j))+
            xlab(paste("average ", namex))+
            ylab(paste("average ", namey))
          #xlim(c(0,1))+
          #ylim(c(0,1))
          plot_coll[[paste0(i, "vs", j)]] = pl
        }
      }
    }
  }
  return(cowplot::plot_grid(plotlist = plot_coll, ncol = 3, label_size = 10, labels = labels))
}

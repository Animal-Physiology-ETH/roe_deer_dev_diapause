# clustering_kmeans
#
#
# selected_mat =



#' @title clustering_kmeans
#' @description
#' get kmeans clusters for genes
#' @param exprmat matrix with data to base clustering on (gene names should be row names)
#' @param selection selection of genes to cluster
#' @param embryo_i matrix with info on embryos
#' @param embryo_ID_column column name for identifiers for embryos
#' @param gene_ID_column column name for identifiers for genes
#'
#' @param seed seed for kmeans
#' @param rowkm number of clusters in rows
#' @param row_km_repeats number of repeats for rows for kmeans
#' @param columnkm number of clusters in columns
#' @param column_km_repeats number of repeats for columns for kmeans
#' @param save save a list object with clustering
#' @param output output file for save
#' @import ComplexHeatmap
#' @export

clustering_kmeans = function(exprmat = exprmat_TMM_minmax,
                        selection = markergenes$GeneName,
                        embryo_i = embryo_info,
                        embryo_ID_column = "E_ID",
                        gene_ID_column = "GeneName",
                        seed = 123,
                        rowkm = 6,
                        columnkm = 4,
                        rowkm_rep = 800,
                        columnkm_rep = 800,
                        save =F,
                        output = file.path("data", "clusterlist.rda")){
  clustering = list()
  exprdf = as.data.frame(t(exprmat))
  selection = selection[selection %in% colnames(exprdf)]
  selected_df = as.matrix(t(exprdf %>%
                    select(all_of(selection))%>%
                    filter(rownames(exprdf) %in% embryo_i$E_ID)
  ))
  set.seed(seed)
  ht = ComplexHeatmap::Heatmap(selected_df,
                               row_km = rowkm,
                               row_km_repeats = rowkm_rep,
                               show_row_dend = FALSE,
                               cluster_rows = F,
                               cluster_row_slices = F,

                               column_km = columnkm,
                               column_km_repeats = columnkm_rep,
                               show_column_dend = FALSE,
                               show_column_names = FALSE,
                               name = "scaled\nTMM\ncounts"
                               )



  set.seed(seed)
  ht = draw(ht)

  gene_cluster = setNames(data.frame(matrix(ncol = 2, nrow = 0)), c(gene_ID_column, "ClusterNo"))

  for (i in 1:length(row_order(ht))) {
    gene_cluster = rbind(gene_cluster,
                         data.frame(GeneName = t(t(row.names(selected_df[row_order(ht)[[i]],]))),
                                    ClusterNo = as.character(rep(LETTERS[i], length(row_order(ht)[[i]])))))
  }
  clustering[["Genes"]] = gene_cluster

  embryo_cluster = setNames(data.frame(matrix(ncol = 2, nrow = 0)), c(embryo_ID_column, "ClusterNo"))

  for (i in 1:length(column_order(ht))) {
    embryo_cluster = rbind(embryo_cluster,
                           data.frame(EmbryoName = t(t(row.names(t(selected_df)[column_order(ht)[[i]],]))),
                                      ClusterNo = as.character(rep(i, length(column_order(ht)[[i]])))))
  }
  embryo_cluster <- embryo_cluster %>%
    mutate(ClusterNo = recode(ClusterNo, "1"="3", "2" = "4", "3" = "1", "4"="2"))
  colnames(embryo_cluster) = c(embryo_ID_column, "ClusterNo")
  clustering[["Embryos"]] = embryo_cluster
  if (save == T) {
    save(clustering, file = output)
  }
  #expr_dfround = round(expr_df, 10)
  return(clustering)
}

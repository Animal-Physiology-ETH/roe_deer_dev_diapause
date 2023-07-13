#' @title std_error
#' @description
#' standard error calculation
#' @export

std_error = function(x) sd(x)/sqrt(length(x))


#' @title star_p
#' @description
#' converts numeric p values to stars
#' @export

star_p = function(x){
  if (x>=0.05) {
    return("ns")
  }
  if (x<0.001) {
    return("***")
  }
  if (x<0.01) {
    return("**")
  }
  if (x<0.05) {
    return("*")
  }
}


#' @title p_value_df
#' @description return df with p value
#' @param df data frame to create p_value data frame from
#' @param variable variable to analyze
#' @param grouping  grouping variable
#' @param log_scale should data be log_tranformed (T/F)
#' @export
#'


p_value_df = function(df,
                      variable,
                      grouping,
                      log_scale = TRUE) {
  n_groups =  length(unique(df[[grouping]]))
  p_value_DF = data.frame(group1 = as.character(seq(1, n_groups-1)),
                          group2 = as.character(seq(2, n_groups)))
  p_value_DF$comparison = paste0(p_value_DF$group1, "v", p_value_DF$group2)
  results =  signif(pairwise.wilcox.test(df[[variable]], df[[grouping]], paired = FALSE,
                                         p.adjust.method = "BH")[["p.value"]], digits = 3)
  p =c()
  y.pos = c()
  for (i in seq(1, n_groups-1)) {
    p = c(p, results[i,i])
    y.pos = c(y.pos, ifelse(log_scale == TRUE,
                            1.2*log10(max(max(df[which(df[[grouping]] == i), variable]),
                                          max(df[which(df[[grouping]] == i+1), variable]))),
                            1.2*(max(max(df[which(df[[grouping]] == i), variable]),
                                     max(df[which(df[[grouping]] == i+1), variable])))))
  }
  p_value_DF$p = p
  p_value_DF$y.pos = y.pos
  p_value_DF$stars<- lapply(p_value_DF$p, star_p)
  p_value_DF$y.pos_remns = ifelse(p_value_DF$stars == "ns", NA, p_value_DF$y.pos)
  return(p_value_DF)
}

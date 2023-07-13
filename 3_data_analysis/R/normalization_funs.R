##### vector normalizing functions #####

#' @title minmax
#' @description performs a min-max standardization of x
#' @param x numeric vector, cannot be all 0
minmax = function(x){(x- min(x)) /(max(x)-min(x))}

#' @title zscore
#' @description performs a z-score standardization of x
#' @param x numeric vector, cannot be all 0
zscore = function(x){(x - mean(x)) / sd(x)}


##### dataframe normalizing functions #####

#' @title scale_df
#' @description performs a min-max standardization of a dataframe, accross one direction
#' @param df a numeric dataframe
#' @param direction direction in which the function will be applied. Either rows or cols, default = cols.
#' @param fun function to be applied
#' @export

scale_df = function(df, fun = minmax, direction = "cols"){
  if (direction == "rows") {
    df = t(df)
  }
  df = df[,colSums(df)>0] #this is to remove all 0 rows and prevent NA
  df = as.data.frame(apply(df, 2, fun))
  if (direction == "rows") {
    df = t(df)
  }
  return(df)
}

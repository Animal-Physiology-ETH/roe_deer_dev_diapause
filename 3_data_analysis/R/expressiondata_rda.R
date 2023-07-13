#  Create a list objects with all flavours of data
#' @title expression_data_rda
#' @description create an Rds object or a list with all types of data matrices.
#' @param data_location folder path containing expression data, default is file.path("data", "raw_tables")
#' @param prefix file name prior to input type, default is "RoeDeerEmbryo_counts_"
#' @param inputtypes vector with types of input (must be in filename before suffix), default is c("expected", "TPM", "FPKM")
#' @param suffix suffix of input file, default is ".txt"
#' @param scaled_counts also create scaled list for each normalization, default = F
#' @param scaling_methods which scaling function to use, default = c("minmax", "zscore")
#' @param outdir directory for output
#' @param outfilename file name for output (without suffix)
#' @param save_list save the the a list as Rds file? (T/F)
#' @param return_list return the generated list
#' @param rem_outlier remove one embryo with low quality reads? (T/F)
#' @returns list object with differently normalized and/or scaled expression data
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @export

expression_data_rda = function(data_location = file.path("extdata", "raw_tables"),
                                prefix = "RoeDeerEmbryo_counts_",
                                inputtypes = c("expected", "TPM", "FPKM", "TMM"),
                                suffix =".txt",
                                save_list = F,
                                scaled_counts = F,
                                scaling_methods = c("minmax", "zscore"),
                                outdir = "data",
                                outfilename = "exprmat_",
                                return_list = F,
                                rem_outlier = T){
  expression_ls = list()
  expression_env = new.env()
  for (i in 1:(length(inputtypes))) {
    if (inputtypes[i]=="TMM") {
      filename = paste0(prefix, "expected", suffix)
      read_counts = as.matrix(read.delim(file.path(data_location, filename), row.names="gene_name"))
      if (rem_outlier == T) {
        read_counts=read_counts[,colnames(read_counts)!="E1371_E1"]
      }
      # EdgeR TMM calculation
      y = DGEList(counts=read_counts)
      y =  calcNormFactors(y)
      norm_read_counts =  cpm(y,  normalized.lib.sizes  = TRUE)
      varname =  paste0(outfilename, "TMM")
      expression_ls[[varname]] = norm_read_counts
      expression_env[[varname]] = read_counts
      outfile =  file.path(outdir, paste0(varname, ".rda"))
      save(list=c(varname), envir=expression_env, file = outfile)
    } else{
      filename = paste0(prefix, inputtypes[i], suffix)
      read_counts = as.matrix(read.delim(file.path(data_location, filename), row.names="gene_name"))
      if (rem_outlier == T) {
        read_counts=read_counts[,colnames(read_counts)!="E1371_E1"]
      }
      varname = paste0(outfilename, inputtypes[i])
      expression_ls[[varname]] = read_counts
      expression_env[[varname]] = read_counts
      outfile =  file.path(outdir, paste0(varname, ".rda"))
      save(list=c(varname), envir=expression_env, file = outfile)
    }
  }
  if (scaled_counts == T) {
    for (i in 1:length(expression_ls)) {
      for (method in scaling_methods) {
        data_i = expression_ls[[i]]
        data_i = data_i[order(rownames(data_i)),]
        data_scaled = scale_df(data_i, method, direction = "rows")
        df_name = paste0(names(expression_ls[i]),"_", method)
        expression_ls[[df_name]] = data_scaled
        expression_env[[df_name]] = data_scaled
        outfile =  file.path(outdir, paste0(df_name, ".rda"))
        save(list=c(df_name), envir=expression_env, file = outfile)
      }
      }
  }
  if (save_list == T) {
    save(expression_ls, file = file.path(outdir, paste0(outfilename, "list", ".Rds")))
  }
  if (return_list == T) {
    return(expression_ls)
  }
}



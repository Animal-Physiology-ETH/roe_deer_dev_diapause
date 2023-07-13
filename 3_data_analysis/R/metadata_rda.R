#  Create a list objects with all flavours of data
#' @title metadata_rda
#' @description creates an rda objects of all metadata
#' @param data_location folder path containing metadata files, default is file.path("data", "raw_tables", "metadata")
#' @param suffix suffix of input file, default is ".xlsx"
#' @param returnlist return the list? (T/F)
#' @param outdir directory for output
#' @returns list object with differently normalized and/or scaled expression data
#' @importFrom readxl read_xlsx
#' @export

metadata_rda = function(data_location = file.path("extdata", "raw_tables", "metadata"),
                        suffix = ".xlsx",
                        returnlist = F,
                        outdir = "data"
                        ){
  filels = list.files(data_location)
  var_env = new.env()
  var_list = list()
  for (file in filels) {
      base_name = gsub(suffix, "", file)
      df = read_xlsx(file.path(data_location, file))
      var_list[[paste(base_name)]] = df
      var_env[[paste(base_name)]] = df
      save(list=c(base_name), envir=var_env, file = file.path(outdir, paste0(gsub(suffix, "", file), ".rda")))
  }
  if (returnlist == T) {
    return(var_list)
  }
}

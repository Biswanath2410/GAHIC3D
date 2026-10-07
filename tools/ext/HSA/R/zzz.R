#' @useDynLib "hsa", .registration = TRUE, .fixes = "C_"
.onUnload <- function(libpath) {
  library.dynam.unload("hsa", libpath)
}

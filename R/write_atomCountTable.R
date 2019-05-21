
#' Write default atomCountTable to current working directory
#'
#' @param path file path to write to
#'
#' @return TRUE is sucsssful
#' @export
write_atomCountTable <- function(path = getwd(), overwrite = F){
  base::file.copy(system.file('defaultResidueAtoms.txt', package = 'peptideUtils', mustWork = T),
                  to = path,
                  overwrite = overwrite,
                  recursive = F)
}

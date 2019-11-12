
#' Write default atomCountTable to current working directory
#'
#' @param path file path to write to
#'
#' @return TRUE is successful
#' @export
write_atomCountTable <- function(path = getwd(), overwrite = F){
  base::file.copy(system.file('defaultResidueAtoms.txt', package = 'peptideUtils', mustWork = T),
                  to = path,
                  overwrite = overwrite,
                  recursive = F)
}

#' Read a list of IDs, and protein names from a fasta file.
#' 
#' @param fname Path to fasta file to get data from.
#' @param id_line_regex Regex to find fasta entry header lines
#' @param capture_regex Regex with 2 capturing groups.
#' The first will be put in the \code{id} column, the second will be put in the \code{protein}.
#' @param .colnames Character vector of length 2 containing the column names in the output.
#' 
#' @examples
#' read_ids()
#' 
#' @export
read_ids <- function(fname = system.file('extdata/Human_uniprot-reviewed_20171020.fasta',
	package = 'peptideUtils', mustWork = T),
	id_line_regex = '^>(sp|tr)',
	capture_regex = '^>[sptr]{2}\\|(\\w+)\\|(\\w+)_\\w+',
	.colnames = c('id', 'protein'))
{
	if(length(.colnames) != 2){
		stop('Length of .colnames must be 2!')
	}

	fasta <- scan(fname, what = character(), sep = '\n', quiet = T)
	fasta <- fasta[grepl(id_line_regex, fasta)]
	matches <- regmatches(fasta, regexec(capture_regex, fasta))
	d <- as.data.frame(do.call(rbind, matches), stringsAsFactors = F)
	names(d) <- c('V1', .colnames)
	return(d[,.colnames])
}


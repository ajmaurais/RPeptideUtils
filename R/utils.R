
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
read_ids <- function(
  fname = system.file(
    'extdata/Human_uniprot-reviewed_20171020.fasta',
    package = 'peptideUtils', mustWork = T
  ),
  id_line_regex = '^>(sp|tr)',
  capture_regex = '^>[sptr]{2}\\|(\\w+)\\|(\\w+)_\\w+',
  .colnames = c('id', 'protein')
) {
  if(length(.colnames) != 2){
    top('Length of .colnames must be 2!')
  }

  fasta <- scan(fname, what = character(), sep = '\n', quiet = T)
  fasta <- fasta[grepl(id_line_regex, fasta)]
  matches <- regmatches(fasta, regexec(capture_regex, fasta))
  d <- as.data.frame(do.call(rbind, matches), stringsAsFactors = F)
  names(d) <- c('V1', .colnames)
  d[,.colnames]
}

#' Download protein sequences from the UniProt API for a list of accessions.
#'
#' @param accessions A vector of UniProt accessions to download.
#' @param output_file An optional file path to a fasta file to write sequences to.
#' @param use_cached If TRUE, read and return the matching ids in \code{output_file} if the file exists.
#'        Otherwise download the sequences through the UniProt API.
#'        An error is raised if any accession is not found in the cached file.
#' @param batch_size Number of accessions to fetch per API request (default 100)
#'
#' @return A data.frame with columns \code{id} and \code{sequence},
#'         matching the format of \code{\link{readFasta}}.
#'
#' @examples
#' fetch_uniprot_sequences(c('P00533'))
#'
#' @export
fetch_uniprot_sequences <- function(
  accessions, output_file=NULL,
  use_cached=FALSE, batch_size=100
) {
  if (!is.character(accessions) || length(accessions) == 0) {
    stop("accessions must be a non-empty character vector")
  }

  accessions <- unique(accessions)

  # Return cached results if available
  if (use_cached && !is.null(output_file) && file.exists(output_file)) {
    cached <- readFasta(output_file)
    missing <- setdiff(accessions, cached$id)
    if (length(missing) > 0) {
      stop("Accessions not found in cached file: ", paste(missing, collapse = ", "))
    }
    return(cached[cached$id %in% accessions, , drop = FALSE])
  }

  # Fetch from UniProt API in batches, writing to a temp file
  base_url <- "https://rest.uniprot.org/uniprotkb/stream"
  tmp <- tempfile(fileext = ".fasta")
  on.exit(unlink(tmp), add = TRUE)

  for (i in seq(1, length(accessions), by = batch_size)) {
    batch <- accessions[i:min(i + batch_size - 1, length(accessions))]
    query <- paste(sprintf("accession:%s", batch), collapse = "+OR+")
    request_url <- sprintf("%s?query=(%s)&format=fasta", base_url, query)

    tryCatch({
      fasta_lines <- .fetch_fasta_lines(request_url)
      write(fasta_lines, file = tmp, append = TRUE)
    }, error = function(e) {
      warning(sprintf("Failed to fetch batch starting at index %d: %s", i, e$message))
    })
  }

  result <- readFasta(tmp)

  if (!is.null(output_file)) {
    file.copy(tmp, output_file, overwrite = TRUE)
  }

  result
}

.fetch_fasta_lines <- function(request_url) {
  con <- url(request_url)
  on.exit(close(con))
  readLines(con, warn = FALSE)
}

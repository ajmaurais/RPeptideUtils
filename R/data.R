#' 
#' pritein id, name lookup
#' 
#' data.frame containing protein uniprot IDs and names. Currently only names for human proteins are included. 
#' 
#' The text representaion of the data frame was generated from the default fasta file included in the package with the folowing shell command:
#'
#' \code{$ cat Human_uniprot-reviewed_20171020.fasta |egrep '^>(tr|sp)'|sed -E 's/>[trsp]{2}\|([A-Z0-9]+)\|([0-9A-Z]+)_HUMAN.*/\1	\2/' > id_lookup.tsv}
#' 
#' @format A data frame with two variables:
#' \describe{
#' 	\item{\code{id}}{uniprot id}
#' 	\item{\code{protein}}{Protein name}
#' }
#' 
#' @examples
#' #add a column containing protein names to a data.frame with a column containing uniprot IDs.
#' data(id_protein_lookup)
#' seqs <- data.frame(id = fastaInfo()[['ids']][1:100], stringsAsFactors = F)
#' seqs <- merge(x = seqs, y = id_protein_lookup, by = 'id')
#' 
'id_protein_lookup'
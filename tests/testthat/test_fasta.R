
context('FastaFile')

ids <- c("A0MZ66", "A6NMY6")

sequences <- c(paste0("MNSSDEEKQLQLITSLKEQAIGEYEDLRAENQKTKEKCDKIRQERDEAVKKLEEFQKISHMVIEEVNFMQNHLEIEKTCRESAEALAT",
                        "KLNKENKTLKRISMLYMAKLGPDVITEEINIDDEDSTTDTDGAAETCVSVQCQKQIKELRDQIVSVQEEKKILAIELENLKSKLVEVI",
                        "EEVNKVKQEKTVLNSEVLEQRKVLEKCNRVSMLAVEEYEEMQVNLELEKDLRKKAESFAQEMFIEQNKLKRQSHLLLQSSIPDQQLLK",
                        "ALDENAKLTQQLEEERIQHQQKVKELEEQLENETLHKEIHNLKQQLELLEEDKKELELKYQNSEEKARNLKHSVDELQKRVNQSENSV",
                        "PPPPPPPPPLPPPPPNPIRSLMSMIRKRSHPSGSGAKKEKATQPETTEEVTDLKRQAVEEMMDRIKKGVHLRPVNQTARPKTKPESSK",
                        "GCESAVDELKGILGTLNKSTSSRSLKSLDPENSETELERILRRRKVTAEADSSSPTGILATSESKSMPVLGSVSSVTKTALNKKTLEA",
                        "EFNSPSPPTPEPGEGPRKLEGCTSSKVTFQPPSSIGCRKKYIDGEKQAEPVVVLDPVSTHEPQTKDQVAEKDPTQHKEDEGEIQPENK",
                        "EDSIENVRETDSSNC", collapse = ''), 
                 paste0("MSTVHEILCKLSLEGDHSTPPSAYGSVKAYTNFDAERDALNIETAIKTKGVDEVTIVNIVTNRDNAQRQDIVFSYQRRTKKELASALK",
                        "SALSGHLETVILGLLKTPAQYDASELKASMKGLGTDEDSLIEIICSRTNQELQEINRVYKEMYKTDLEKDIISDTSGDFRKLMVALAK",
                        "GRRAEDGSVIDYELIDQDAQDLYDAGVKRKGTDVPKWISIMTERSVPHLQKVFDRYKSYSPYDMLESIRKEVKGDLENAFLNLVQRIQ",
                        "NKPLYFADQLYDSMKGKGTRDKVLIRIMVSRSEVDMLKIRSEFKRKYGKSLYYYIQQDTKGDYQKALLYLCGGDD", colapse = ''))


test_that('Correct residues are modified', {
  expect_equal(peptideUtils::getModifiedResidues(c("Q00839", "Q9HCS7", "Q7L014"),
                                                 c("APQC*LGK", "FADMEC*K", "GAEIIVC*TPGR"),
                                                 system.file('extdata/Human_uniprot-reviewed_20171020.fasta',
                                                  			 package = 'peptideUtils', mustWork = T)),
               c("C562", "C676", "C501"))
})

test_that('Correct protein sequences', {
  expect_equal(peptideUtils::getSequences(ids,
                                          system.file('extdata/Human_uniprot-reviewed_20171020.fasta',
                                          			  package = 'peptideUtils', mustWork = T)),
  			   sequences)
})

d <- peptideUtils::fastaInfo(system.file('extdata/Human_uniprot-reviewed_20171020.fasta',
                                         package = 'peptideUtils', mustWork = T))

test_that('Fasta metadata is correct.', {
	expect_equal(d$seq_count, length(d$ids))
})

test_that('Read fasta.', {
	expect_error(peptideUtils::readFasta(system.file('extdata/Human_uniprot-reviewed_20171020.fasta',
                                         package = 'peptideUtils', mustWork = T), -1)
	)
})

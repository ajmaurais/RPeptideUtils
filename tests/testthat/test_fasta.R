
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
  expect_equal(RPeptideUtils::getModifiedResidues(c("Q00839", "Q9HCS7", "Q7L014"),
                                                 c("APQC*LGK", "FADMEC*K", "GAEIIVC*TPGR"),
                                                 system.file('extdata/Human_uniprot-reviewed_20171020.fasta',
                                                  			 package = 'RPeptideUtils', mustWork = T)),
               c("C562", "C676", "C501"))
})

test_that('Correct protein sequences', {
  expect_equal(RPeptideUtils::getSequences(ids,
                                          system.file('extdata/Human_uniprot-reviewed_20171020.fasta',
                                          			  package = 'RPeptideUtils', mustWork = T)),
  			   sequences)
})

d <- RPeptideUtils::fastaInfo(system.file('extdata/Human_uniprot-reviewed_20171020.fasta',
                                         package = 'RPeptideUtils', mustWork = T))

test_that('Fasta metadata is correct.', {
	expect_equal(d$seq_count, length(d$ids))
})

test_that('Read fasta.', {
	expect_error(RPeptideUtils::readFasta(system.file('extdata/Human_uniprot-reviewed_20171020.fasta',
                                         package = 'RPeptideUtils', mustWork = T), -1)
	)
})

test_that('nAfter gives corect number of residues', {
    expect_equal(nchar(RPeptideUtils::nAfter(c('KLNKENKTL', 'SALSGHLE'), sequences, 1)), rep(1, 2))
    expect_equal(nchar(RPeptideUtils::nAfter(c('KLNKENKTL', 'SALSGHLE'), sequences, 2)), rep(2, 2))
    expect_equal(nchar(RPeptideUtils::nAfter(c('KLNKENKTL', 'SALSGHLE'), sequences, 3)), rep(3, 2))
})

test_that('nBefore gives corect number of residues', {
    expect_equal(nchar(RPeptideUtils::nBefore(c('KLNKENKTL', 'SALSGHLE'), sequences, 1)), rep(1, 2))
    expect_equal(nchar(RPeptideUtils::nBefore(c('KLNKENKTL', 'SALSGHLE'), sequences, 2)), rep(2, 2))
    expect_equal(nchar(RPeptideUtils::nBefore(c('KLNKENKTL', 'SALSGHLE'), sequences, 3)), rep(3, 2))
})

test_that('nAfter gives corect residues', {
    expect_equal(RPeptideUtils::nAfter(c('KLNKENKTL', 'SALSGHLE'), sequences, 1), c('K', 'T'))
    expect_equal(RPeptideUtils::nAfter(c('KLNKENKTL', 'SALSGHLE'), sequences, 2), c('KR', 'TV'))
    expect_equal(RPeptideUtils::nAfter(c('KLNKENKTL', 'SALSGHLE'), sequences, 3), c('KRI', 'TVI'))
})

test_that('nBefore gives corect number of residues', {
    expect_equal(RPeptideUtils::nBefore(c('KLNKENKTL', 'SALSGHLE'), sequences, 1), c('T', 'K'))
    expect_equal(RPeptideUtils::nBefore(c('KLNKENKTL', 'SALSGHLE'), sequences, 2), c('AT', 'LK'))
    expect_equal(RPeptideUtils::nBefore(c('KLNKENKTL', 'SALSGHLE'), sequences, 3), c('LAT', 'ALK'))
})

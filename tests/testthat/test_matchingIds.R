
context('matchingProteins')

peptides = c('NHFEPGVYVCAKCGY', 'FFADVWRICTNNTNCTVIND', 'SPQGLKPVVSTEAPPIIFA', 'PEPTI')

ids = list('NHFEPGVYVCAKCGY' = 'Q9NZV6',
           'FFADVWRICTNNTNCTVIND' = 'P54851',
           'SPQGLKPVVSTEAPPIIFA' = 'O14548',
           'PEPTI' = c('Q96RW7', 'P10721', 'Q8IW35', 'Q9HCK4', 'Q9Y6N7', 'Q6ZNJ1'))

multiMatches = c('A' = 20174, 'C' = 19657, 'D' = 20072, 'E' = 20151, 'F' = 20016, 'G' = 20189,
                 'H' = 19776, 'I' = 20033, 'K' = 20114, 'L' = 20167, 'M' = 20213, 'N' = 19957,
                 'O' = 0, 'P' = 20159, 'Q' = 20141, 'R' = 20160, 'S' = 20210, 'T' = 20167,
                 'U' = 26, 'V' = 20173, 'W' = 18924, 'X' = 1, 'Y' = 19802,
                 'AA' = 16580, 'CC' = 5014, 'FART' = 64, 'PEE' = 3967, 'ASS' = 4977)


sortList <- function(l) {
    l <- l[sort(names(l))]
    lapply(l, sort)
}

test_that('Results are correct', {
    expect_equal(sortList(RPeptideUtils::matchingProteins(peptides,
                                                          fastaPath = system.file('extdata/Human_uniprot-reviewed_20171020.fasta',
                                                                                  package = 'RPeptideUtils', mustWork = T),
                                                          progressBar = F, nThread = 1)),
                 sortList(ids))
})

test_that('Parallel processing works', {
    expect_equal(sortList(RPeptideUtils::matchingProteins(peptides,
                                                          fastaPath = system.file('extdata/Human_uniprot-reviewed_20171020.fasta',
                                                                                  package = 'RPeptideUtils', mustWork = T),
                                                          progressBar = F, nThread = 2)),
                 sortList(ids))
    expect_equal(sortList(RPeptideUtils::matchingProteins(peptides,
                                            fastaPath = system.file('extdata/Human_uniprot-reviewed_20171020.fasta',
                                                                    package = 'RPeptideUtils', mustWork = T),
                                            progressBar = F)),
                 sortList(ids))
})

test_that("Edge cases don't fail", {
    # A peptide that doesn't exist in the file
    expect_equal(RPeptideUtils::matchingProteins('THISPEPTIDEDSENTEXIST',
                                                 fastaPath = system.file('extdata/Human_uniprot-reviewed_20171020.fasta',
                                                                         package = 'RPeptideUtils', mustWork = T),
                                                 nThread = 1, progressBar = F),
                list('THISPEPTIDEDSENTEXIST' = character(0)))
    # peptides with a lot of matches
    expect_equal(sortList(lapply(RPeptideUtils::matchingProteins(names(multiMatches),
                                                                 fastaPath = system.file('extdata/Human_uniprot-reviewed_20171020.fasta',
                                                                                         package = 'RPeptideUtils', mustWork = T),
                                                                 nThread = 1, progressBar = F), length)),
                 sortList(multiMatches))
})


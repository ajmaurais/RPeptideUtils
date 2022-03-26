
context('getFragmentIons')

fragment_ions = c('b1'='P', 'b2'='PE', 'b3'='PEP', 'b4'='PEPT', 'b5'='PEPTI', 'b6'='PEPTID', 'b7'='PEPTIDE',
                  'y2'='DE', 'y1'='E', 'y3'='IDE', 'y4'='TIDE', 'y5'='PTIDE', 'y6'='EPTIDE', 'y7'='PEPTIDE')

test_that('Correct fragments are calculated', {
    expect_equal(sort(fragment_ions),
                 sort(apply(RPeptideUtils::getFragmentIonIndices('PEPTIDE'), 1, function(x) substr('PEPTIDE', x[1] + 1, x[2] + x[1]))))
    expect_equal(sort(fragment_ions), sort(RPeptideUtils::getFragmentIonSequences('PEPTIDE')))
})



context('nMissed')

test_that('Zero missed cleavages', {
    # Peptides ending in R or K with no internal cleavage sites
    expect_equal(nMissed(c('PEPTIDER', 'PEPTIDEK', 'AAAAAA')), c(0L, 0L, 0L))
})

test_that('One missed cleavage', {
    # K followed by non-P (internal)
    expect_equal(nMissed('PEPKTIDER'), 1L)
    # R followed by non-P (internal)
    expect_equal(nMissed('PEPRTIDER'), 1L)
})

test_that('Multiple missed cleavages', {
    # Two missed cleavages
    expect_equal(nMissed('PEPKTIDERK'), 2L)
    # Three missed cleavages
    expect_equal(nMissed('PEPKTIDERKAK'), 3L)
})

test_that('KP and RP do not count as cleavage sites', {
    # K followed by P should not count
    expect_equal(nMissed('PEPKPIDER'), 0L)
    # R followed by P should not count
    expect_equal(nMissed('PEPRPIDER'), 0L)
    # Mix of cleavable and non-cleavable
    expect_equal(nMissed('PEPKPIDEKR'), 1L)
})

test_that('Vector input works', {
    seqs <- c('PEPTIDER', 'PEPKTIDER', 'PEPKTIDERK', 'PEPKPIDER')
    expect_equal(nMissed(seqs), c(0L, 1L, 2L, 0L))
})

test_that('Custom cleavage pattern works', {
    # Use pattern that only matches K (not R)
    expect_equal(nMissed('PEPKTIDER', cleavagePattern = '(K)(?=[^P])'), 1L)
    expect_equal(nMissed('PEPRTIDER', cleavagePattern = '(K)(?=[^P])'), 0L)
})

test_that('Edge cases', {
    # Single residue
    expect_equal(nMissed('K'), 0L)
    expect_equal(nMissed('R'), 0L)
    # Two residues
    expect_equal(nMissed('KA'), 1L)
    expect_equal(nMissed('KP'), 0L)
    # Empty string
    expect_equal(nMissed(''), 0L)
})

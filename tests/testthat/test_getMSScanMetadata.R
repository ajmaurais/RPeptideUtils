
context('getMSScanMetadata')

mzml_path <- system.file('extdata/20240927_STEL1_Evo3_MW_15N_MS3_dilution_12_1.mzML',
                          package = 'RPeptideUtils', mustWork = TRUE)
mzml_basename <- '20240927_STEL1_Evo3_MW_15N_MS3_dilution_12_1'

test_that('basic scan metadata is correct', {
    result <- RPeptideUtils::getMSScanMetadata(c(1L, 50L, 138L), rep(mzml_path, 3), progressBar = FALSE)

    expect_equal(nrow(result), 3)
    expect_equal(result$scan, c(1L, 50L, 138L))
    expect_equal(result$file, rep(mzml_basename, 3))
    expect_equal(result$path, rep(mzml_path, 3))
    expect_equal(result$ms_level, rep(2L, 3))
})

test_that('retention times are correct', {
    result <- RPeptideUtils::getMSScanMetadata(c(1L, 138L), rep(mzml_path, 2), progressBar = FALSE)

    # Scan 1 RT: 7.267534e-04 minutes = 0.04360520 seconds
    expect_equal(result$rt[1], 0.04360520, tolerance = 1e-6)
    # Scan 138 RT
    expect_equal(result$rt[2], 9.872503, tolerance = 1e-4)
})

test_that('isolation window values are correct', {
    result <- RPeptideUtils::getMSScanMetadata(c(1L, 138L), rep(mzml_path, 2), progressBar = FALSE)

    # Scan 1: target=406.0, offsets=25.5
    expect_equal(result$precursor_mz[1], 406.0)
    expect_equal(result$isolation_window_lower_offset[1], 25.5)
    expect_equal(result$isolation_window_upper_offset[1], 25.5)

    # Scan 138: target=391.9545898, offsets=0.5
    expect_equal(result$precursor_mz[2], 391.9546, tolerance = 1e-4)
    expect_equal(result$isolation_window_lower_offset[2], 0.5)
    expect_equal(result$isolation_window_upper_offset[2], 0.5)
})

test_that('output has expected columns', {
    result <- RPeptideUtils::getMSScanMetadata(1L, mzml_path, progressBar = FALSE)
    expected_cols <- c('scan', 'file', 'path', 'rt', 'ms_level',
                       'precursor_mz', 'isolation_window_lower_offset',
                       'isolation_window_upper_offset')
    expect_equal(colnames(result), expected_cols)
})

test_that('empty input returns empty dataframe', {
    result <- RPeptideUtils::getMSScanMetadata(integer(0), character(0), progressBar = FALSE)
    expect_equal(nrow(result), 0)
    expect_true(is.data.frame(result))
})

test_that('input order is preserved', {
    result <- RPeptideUtils::getMSScanMetadata(c(138L, 1L, 50L), rep(mzml_path, 3), progressBar = FALSE)
    expect_equal(result$scan, c(138L, 1L, 50L))
})

test_that('mismatched input lengths throw error', {
    expect_error(RPeptideUtils::getMSScanMetadata(c(1L, 2L), mzml_path, progressBar = FALSE))
})

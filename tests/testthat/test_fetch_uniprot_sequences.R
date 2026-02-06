
context('fetch_uniprot_sequences')

# ---------------------------------------------------------------------------
# Reference data from the FASTA file included in the package
# ---------------------------------------------------------------------------
fasta_path <- system.file('extdata/Human_uniprot-reviewed_20171020.fasta',
                          package = 'RPeptideUtils', mustWork = TRUE)
ref <- RPeptideUtils::readFasta(fasta_path)
all_lines <- readLines(fasta_path, warn = FALSE)
header_indices <- grep("^>", all_lines)

# Three accessions used across the tests
test_ids <- c("A0PJW8", "P58062", "Q14493")
test_ref <- ref[match(test_ids, ref$id), ]

# ---------------------------------------------------------------------------
# Helper: extract raw FASTA lines for a set of accessions
# ---------------------------------------------------------------------------
.extract_fasta_entries <- function(ids) {
  out <- character(0)
  for (id in ids) {
    idx <- grep(paste0("\\|", id, "\\|"), all_lines[header_indices])
    if (length(idx) > 0) {
      start <- header_indices[idx]
      end <- if (idx < length(header_indices)) header_indices[idx + 1] - 1 else length(all_lines)
      out <- c(out, all_lines[start:end])
    }
  }
  out
}

# ---------------------------------------------------------------------------
# Mock for .fetch_fasta_lines: parses the request URL to determine which
# accessions were requested and returns only those FASTA lines.
# ---------------------------------------------------------------------------
.mock_fetch <- function(request_url) {
  matches <- regmatches(request_url, gregexpr("accession:(\\w+)", request_url))[[1]]
  requested <- sub("accession:", "", matches)
  .extract_fasta_entries(requested)
}

# ===========================================================================
# Input validation
# ===========================================================================
test_that('errors on non-character input', {
  expect_error(RPeptideUtils::fetch_uniprot_sequences(123),
               "non-empty character vector")
})

test_that('errors on empty character vector', {
  expect_error(RPeptideUtils::fetch_uniprot_sequences(character(0)),
               "non-empty character vector")
})

# ===========================================================================
# Basic fetch
# ===========================================================================
test_that('returns data.frame with id and sequence columns', {
  local_mocked_bindings(.fetch_fasta_lines = .mock_fetch, .package = "RPeptideUtils")
  result <- RPeptideUtils::fetch_uniprot_sequences(test_ids)

  expect_s3_class(result, "data.frame")
  expect_named(result, c("id", "sequence"))
  expect_equal(nrow(result), length(test_ids))
})

test_that('returned sequences match readFasta reference', {
  local_mocked_bindings(.fetch_fasta_lines = .mock_fetch, .package = "RPeptideUtils")
  result <- RPeptideUtils::fetch_uniprot_sequences(test_ids)

  for (id in test_ids) {
    expect_equal(result$sequence[result$id == id],
                 test_ref$sequence[test_ref$id == id])
  }
})

# ===========================================================================
# Deduplication
# ===========================================================================
test_that('duplicate accessions are deduplicated', {
  local_mocked_bindings(.fetch_fasta_lines = .mock_fetch, .package = "RPeptideUtils")
  result <- RPeptideUtils::fetch_uniprot_sequences(c("A0PJW8", "A0PJW8", "P58062"))

  expect_equal(nrow(result), 2)
  expect_true(all(c("A0PJW8", "P58062") %in% result$id))
})

# ===========================================================================
# Batching
# ===========================================================================
test_that('batch_size controls number of API requests', {
  call_count <- 0
  local_mocked_bindings(
    .fetch_fasta_lines = function(request_url) {
      call_count <<- call_count + 1
      .mock_fetch(request_url)
    },
    .package = "RPeptideUtils"
  )

  result <- RPeptideUtils::fetch_uniprot_sequences(test_ids, batch_size = 1)

  expect_equal(call_count, length(test_ids))
  expect_equal(nrow(result), length(test_ids))
})

test_that('single batch when batch_size >= n accessions', {
  call_count <- 0
  local_mocked_bindings(
    .fetch_fasta_lines = function(request_url) {
      call_count <<- call_count + 1
      .mock_fetch(request_url)
    },
    .package = "RPeptideUtils"
  )

  result <- RPeptideUtils::fetch_uniprot_sequences(test_ids, batch_size = 100)

  expect_equal(call_count, 1)
  expect_equal(nrow(result), length(test_ids))
})

# ===========================================================================
# output_file
# ===========================================================================
test_that('output_file is written and readable by readFasta', {
  tmp <- tempfile(fileext = ".fasta")
  on.exit(unlink(tmp), add = TRUE)

  local_mocked_bindings(.fetch_fasta_lines = .mock_fetch, .package = "RPeptideUtils")
  result <- RPeptideUtils::fetch_uniprot_sequences(test_ids, output_file = tmp)

  expect_true(file.exists(tmp))
  cached <- RPeptideUtils::readFasta(tmp)
  expect_equal(nrow(cached), nrow(result))
  expect_equal(sort(cached$id), sort(result$id))
})

# ===========================================================================
# use_cached
# ===========================================================================
test_that('use_cached returns cached results without calling API', {
  tmp <- tempfile(fileext = ".fasta")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(.extract_fasta_entries(test_ids), tmp)

  local_mocked_bindings(
    .fetch_fasta_lines = function(...) stop("API should not be called when using cache"),
    .package = "RPeptideUtils"
  )

  result <- RPeptideUtils::fetch_uniprot_sequences(test_ids, output_file = tmp, use_cached = TRUE)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), length(test_ids))
  for (id in test_ids) {
    expect_equal(result$sequence[result$id == id],
                 test_ref$sequence[test_ref$id == id])
  }
})

test_that('use_cached returns only the requested accessions', {
  tmp <- tempfile(fileext = ".fasta")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(.extract_fasta_entries(test_ids), tmp)

  result <- RPeptideUtils::fetch_uniprot_sequences("A0PJW8", output_file = tmp, use_cached = TRUE)

  expect_equal(nrow(result), 1)
  expect_equal(result$id, "A0PJW8")
})

test_that('use_cached errors when accession is missing from cache', {
  tmp <- tempfile(fileext = ".fasta")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(.extract_fasta_entries("A0PJW8"), tmp)

  expect_error(
    RPeptideUtils::fetch_uniprot_sequences(c("A0PJW8", "XXXXXX"),
                                           output_file = tmp, use_cached = TRUE),
    "not found in cached file"
  )
})

test_that('use_cached without output_file falls through to API', {
  local_mocked_bindings(.fetch_fasta_lines = .mock_fetch, .package = "RPeptideUtils")
  result <- RPeptideUtils::fetch_uniprot_sequences(test_ids, use_cached = TRUE)

  expect_equal(nrow(result), length(test_ids))
})

test_that('use_cached with nonexistent output_file falls through to API', {
  tmp <- tempfile(fileext = ".fasta")
  on.exit(unlink(tmp), add = TRUE)

  local_mocked_bindings(.fetch_fasta_lines = .mock_fetch, .package = "RPeptideUtils")
  result <- RPeptideUtils::fetch_uniprot_sequences(test_ids, output_file = tmp, use_cached = TRUE)

  expect_equal(nrow(result), length(test_ids))
  expect_true(file.exists(tmp))
})

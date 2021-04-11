
context("oneLetterToThree")

one_letter_sequences = c("ACLLPETVNMEEYPYDAEY", "ALC*AEFK", "AQUPIVER", "CTGGEVGATSALAPK", "IVSNASCTTNCLAPLAK",
                        "LDADIHTNTCR", "LLYVSCNPR", "LPACVVDCGTGYTK", "MIVECVMNNATCTR", "SEGLPSECR", "TPCGEGSK")

three_letter_sequences = c("AlaCysLeuLeuProGluThrValAsnMetGluGluTyrProTyrAspAlaGluTyr", "AlaLeuCys*AlaGluPheLys",
                           "AlaGlnSecProIleValGluArg", "CysThrGlyGlyGluValGlyAlaThrSerAlaLeuAlaProLys",
                           "IleValSerAsnAlaSerCysThrThrAsnCysLeuAlaProLeuAlaLys", "LeuAspAlaAspIleHisThrAsnThrCysArg",
                           "LeuLeuTyrValSerCysAsnProArg", "LeuProAlaCysValValAspCysGlyThrGlyTyrThrLys",
                           "MetIleValGluCysValMetAsnAsnAlaThrCysThrArg", "SerGluGlyLeuProSerGluCysArg", "ThrProCysGlyGluGlySerLys")

test_that("one letter conversions work", {
  expect_equal(oneLetterToThree(one_letter_sequences), three_letter_sequences)
  expect_error(oneLetterToThree('*FKNLR'))
  expect_error(oneLetterToThree('FKNBLR'))
  expect_error(oneLetterToThree('PheLysAsnAssLeuArg'))
})

test_that("deliminators are handeled", {
  expect_equal(oneLetterToThree('ALC*AEFK', sep_out = '-'), "Ala-Leu-Cys*-Ala-Glu-Phe-Lys")
  expect_equal(oneLetterToThree('LLYVSCNPR', sep_out = '-'), "Leu-Leu-Tyr-Val-Ser-Cys-Asn-Pro-Arg")
  expect_equal(oneLetterToThree('A-L-C*-A-E-F-K', sep_in = '-', sep_out = '-'), "Ala-Leu-Cys*-Ala-Glu-Phe-Lys")
  expect_equal(oneLetterToThree('L-L-Y-V-S-C-N-P-R', sep_in = '-', sep_out = '-'), "Leu-Leu-Tyr-Val-Ser-Cys-Asn-Pro-Arg")
  expect_equal(oneLetterToThree('ALC*AEFK', sep_out = '-', n_term_out = 'H', c_term_out = 'OH'), "H-Ala-Leu-Cys*-Ala-Glu-Phe-Lys-OH")
  expect_equal(oneLetterToThree('LLYVSCNPR', sep_out = '-', n_term_out = 'H', c_term_out = 'OH'), "H-Leu-Leu-Tyr-Val-Ser-Cys-Asn-Pro-Arg-OH")
})

context("threeLetterToOne")

test_that("three letter conversions work", {
  expect_equal(threeLetterToOne(three_letter_sequences), one_letter_sequences)
  expect_error(threeLetterToOne('*PheLysAsnLeuArg'))
  expect_error(threeLetterToOne('PheLysAsnAssLeuArg'))
  expect_error(threeLetterToOne('FKNBLR'))
})

test_that("deliminators are handeled", {
  expect_equal(threeLetterToOne("Ala-Leu-Cys*-Ala-Glu-Phe-Lys", sep_in = '-'), 'ALC*AEFK')
  expect_equal(threeLetterToOne("Leu-Leu-Tyr-Val-Ser-Cys-Asn-Pro-Arg", sep_in = '-'), 'LLYVSCNPR')
  expect_equal(threeLetterToOne("Ala-Leu-Cys*-Ala-Glu-Phe-Lys", sep_in = '-', sep_out = '-'), 'A-L-C*-A-E-F-K')
  expect_equal(threeLetterToOne("Leu-Leu-Tyr-Val-Ser-Cys-Asn-Pro-Arg", sep_in = '-', sep_out = '-'), 'L-L-Y-V-S-C-N-P-R')
})

context("Functons interconvert")

test_that("functions can convert back and forth", {
  expect_equal(threeLetterToOne(oneLetterToThree(one_letter_sequences)), one_letter_sequences)
  expect_equal(oneLetterToThree(threeLetterToOne(three_letter_sequences)), three_letter_sequences)
})

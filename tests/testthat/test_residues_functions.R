
context("calcMass")

standards <- data.frame(sequence = c("ACLLPETVNMEEYPYDAEY", "ALCAEFK", "AQCPIVER", "CTGGEVGATSALAPK", "IVSNASCTTNCLAPLAK",
                                     "LDADIHTNTCR", "LLYVSCNPR", "LPACVVDCGTGYTK", "MIVECVMNNATCTR", "SEGLPSECR", "TPCGEGSK"),
                        formula = c("C102H147N21O36S2", "C37H59N9O11S", "C40H69N13O13S", "C58H99N17O22S",
                                    "C75H130N22O26S2", "C52H86N18O20S", "C49H80N14O14S", "C65H105N17O22S2",
                                    "C66H115N21O23S4", "C40H67N13O17S", "C32H54N10O14S"),
                        exact_mass = c(2305.9764, 837.4057, 971.4861, 1417.6874, 1818.8972, 1314.5989, 1120.5702, 
                                       1539.7065, 1697.7361, 1033.4501, 834.3543),
                        molecular_weight = c(2307.529, 837.991, 972.130, 1418.587, 1820.113, 1315.426, 1121.323,
                                             1540.772, 1699.010, 1034.110, 834.900),
                        stringsAsFactors = F)

test_that("Correct masses are calculated",{
  
  #check exact mass against standards
  expect_equal(peptideUtils::calcMass(standards$sequence,
                                      residueAtoms = system.file('defaultResidueAtoms.txt', package = 'peptideUtils', mustWork = T)),
               standards$exact_mass, 
               tolerance = 1e-3)
  
  #check avg mass against standards
  expect_equal(peptideUtils::calcMass(standards$sequence,
                                      residueAtoms = system.file('defaultResidueAtoms.txt', package = 'peptideUtils', mustWork = T),
                                      monoMass = F),
               standards$molecular_weight, 
               tolerance = 1e-3)
})

context('calcFormula')

test_that("Correct molecular formulas are calculated",{
  
  #check exact mass against standards
  expect_equal(peptideUtils::calcFormula(standards$sequence,
                                         residueAtoms = system.file('defaultResidueAtoms.txt', package = 'peptideUtils', mustWork = T)),
               standards$formula)
})


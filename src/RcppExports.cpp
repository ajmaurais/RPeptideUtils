// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// getSquences
Rcpp::CharacterVector getSquences(const Rcpp::CharacterVector& ids, std::string fastaPath);
RcppExport SEXP _peptideUtils_getSquences(SEXP idsSEXP, SEXP fastaPathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ids(idsSEXP);
    Rcpp::traits::input_parameter< std::string >::type fastaPath(fastaPathSEXP);
    rcpp_result_gen = Rcpp::wrap(getSquences(ids, fastaPath));
    return rcpp_result_gen;
END_RCPP
}
// getModifiedResidues
Rcpp::CharacterVector getModifiedResidues(const Rcpp::CharacterVector& ids, const Rcpp::CharacterVector& peptideSeq, std::string fastaPath, std::string modSep);
RcppExport SEXP _peptideUtils_getModifiedResidues(SEXP idsSEXP, SEXP peptideSeqSEXP, SEXP fastaPathSEXP, SEXP modSepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ids(idsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type peptideSeq(peptideSeqSEXP);
    Rcpp::traits::input_parameter< std::string >::type fastaPath(fastaPathSEXP);
    Rcpp::traits::input_parameter< std::string >::type modSep(modSepSEXP);
    rcpp_result_gen = Rcpp::wrap(getModifiedResidues(ids, peptideSeq, fastaPath, modSep));
    return rcpp_result_gen;
END_RCPP
}
// combineMods
std::string combineMods(const Rcpp::CharacterVector& mods, char sep);
RcppExport SEXP _peptideUtils_combineMods(SEXP modsSEXP, SEXP sepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type mods(modsSEXP);
    Rcpp::traits::input_parameter< char >::type sep(sepSEXP);
    rcpp_result_gen = Rcpp::wrap(combineMods(mods, sep));
    return rcpp_result_gen;
END_RCPP
}
// calcMass
Rcpp::NumericVector calcMass(const Rcpp::StringVector& sequences, bool monoMass, std::string residueAtoms, std::string atomMasses);
RcppExport SEXP _peptideUtils_calcMass(SEXP sequencesSEXP, SEXP monoMassSEXP, SEXP residueAtomsSEXP, SEXP atomMassesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type sequences(sequencesSEXP);
    Rcpp::traits::input_parameter< bool >::type monoMass(monoMassSEXP);
    Rcpp::traits::input_parameter< std::string >::type residueAtoms(residueAtomsSEXP);
    Rcpp::traits::input_parameter< std::string >::type atomMasses(atomMassesSEXP);
    rcpp_result_gen = Rcpp::wrap(calcMass(sequences, monoMass, residueAtoms, atomMasses));
    return rcpp_result_gen;
END_RCPP
}
// calcFormula
Rcpp::StringVector calcFormula(const Rcpp::StringVector& sequences, bool subscripts, std::string residueAtoms);
RcppExport SEXP _peptideUtils_calcFormula(SEXP sequencesSEXP, SEXP subscriptsSEXP, SEXP residueAtomsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type sequences(sequencesSEXP);
    Rcpp::traits::input_parameter< bool >::type subscripts(subscriptsSEXP);
    Rcpp::traits::input_parameter< std::string >::type residueAtoms(residueAtomsSEXP);
    rcpp_result_gen = Rcpp::wrap(calcFormula(sequences, subscripts, residueAtoms));
    return rcpp_result_gen;
END_RCPP
}
// oneLetterToThree
Rcpp::StringVector oneLetterToThree(Rcpp::StringVector sequences, std::string sep_in, std::string sep_out, std::string n_term_out, std::string c_term_out);
RcppExport SEXP _peptideUtils_oneLetterToThree(SEXP sequencesSEXP, SEXP sep_inSEXP, SEXP sep_outSEXP, SEXP n_term_outSEXP, SEXP c_term_outSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type sequences(sequencesSEXP);
    Rcpp::traits::input_parameter< std::string >::type sep_in(sep_inSEXP);
    Rcpp::traits::input_parameter< std::string >::type sep_out(sep_outSEXP);
    Rcpp::traits::input_parameter< std::string >::type n_term_out(n_term_outSEXP);
    Rcpp::traits::input_parameter< std::string >::type c_term_out(c_term_outSEXP);
    rcpp_result_gen = Rcpp::wrap(oneLetterToThree(sequences, sep_in, sep_out, n_term_out, c_term_out));
    return rcpp_result_gen;
END_RCPP
}
// threeLetterToOne
Rcpp::StringVector threeLetterToOne(Rcpp::StringVector sequences, std::string sep_in, std::string sep_out, std::string n_term_out, std::string c_term_out);
RcppExport SEXP _peptideUtils_threeLetterToOne(SEXP sequencesSEXP, SEXP sep_inSEXP, SEXP sep_outSEXP, SEXP n_term_outSEXP, SEXP c_term_outSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type sequences(sequencesSEXP);
    Rcpp::traits::input_parameter< std::string >::type sep_in(sep_inSEXP);
    Rcpp::traits::input_parameter< std::string >::type sep_out(sep_outSEXP);
    Rcpp::traits::input_parameter< std::string >::type n_term_out(n_term_outSEXP);
    Rcpp::traits::input_parameter< std::string >::type c_term_out(c_term_outSEXP);
    rcpp_result_gen = Rcpp::wrap(threeLetterToOne(sequences, sep_in, sep_out, n_term_out, c_term_out));
    return rcpp_result_gen;
END_RCPP
}
// digest
Rcpp::List digest(Rcpp::CharacterVector sequences, Rcpp::CharacterVector ids, unsigned nMissedCleavages, std::string cleavagePattern, bool mz_filter, std::string residueAtoms, std::string atomMasses, double minMz, double maxMz, int minCharge, int maxCharge, size_t minLen, size_t maxLen);
RcppExport SEXP _peptideUtils_digest(SEXP sequencesSEXP, SEXP idsSEXP, SEXP nMissedCleavagesSEXP, SEXP cleavagePatternSEXP, SEXP mz_filterSEXP, SEXP residueAtomsSEXP, SEXP atomMassesSEXP, SEXP minMzSEXP, SEXP maxMzSEXP, SEXP minChargeSEXP, SEXP maxChargeSEXP, SEXP minLenSEXP, SEXP maxLenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type sequences(sequencesSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type ids(idsSEXP);
    Rcpp::traits::input_parameter< unsigned >::type nMissedCleavages(nMissedCleavagesSEXP);
    Rcpp::traits::input_parameter< std::string >::type cleavagePattern(cleavagePatternSEXP);
    Rcpp::traits::input_parameter< bool >::type mz_filter(mz_filterSEXP);
    Rcpp::traits::input_parameter< std::string >::type residueAtoms(residueAtomsSEXP);
    Rcpp::traits::input_parameter< std::string >::type atomMasses(atomMassesSEXP);
    Rcpp::traits::input_parameter< double >::type minMz(minMzSEXP);
    Rcpp::traits::input_parameter< double >::type maxMz(maxMzSEXP);
    Rcpp::traits::input_parameter< int >::type minCharge(minChargeSEXP);
    Rcpp::traits::input_parameter< int >::type maxCharge(maxChargeSEXP);
    Rcpp::traits::input_parameter< size_t >::type minLen(minLenSEXP);
    Rcpp::traits::input_parameter< size_t >::type maxLen(maxLenSEXP);
    rcpp_result_gen = Rcpp::wrap(digest(sequences, ids, nMissedCleavages, cleavagePattern, mz_filter, residueAtoms, atomMasses, minMz, maxMz, minCharge, maxCharge, minLen, maxLen));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_peptideUtils_getSquences", (DL_FUNC) &_peptideUtils_getSquences, 2},
    {"_peptideUtils_getModifiedResidues", (DL_FUNC) &_peptideUtils_getModifiedResidues, 4},
    {"_peptideUtils_combineMods", (DL_FUNC) &_peptideUtils_combineMods, 2},
    {"_peptideUtils_calcMass", (DL_FUNC) &_peptideUtils_calcMass, 4},
    {"_peptideUtils_calcFormula", (DL_FUNC) &_peptideUtils_calcFormula, 3},
    {"_peptideUtils_oneLetterToThree", (DL_FUNC) &_peptideUtils_oneLetterToThree, 5},
    {"_peptideUtils_threeLetterToOne", (DL_FUNC) &_peptideUtils_threeLetterToOne, 5},
    {"_peptideUtils_digest", (DL_FUNC) &_peptideUtils_digest, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_peptideUtils(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

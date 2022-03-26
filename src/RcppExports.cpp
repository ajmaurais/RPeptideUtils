// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// getFragmentIonIndices
Rcpp::IntegerMatrix getFragmentIonIndices(std::string seq);
RcppExport SEXP _RPeptideUtils_getFragmentIonIndices(SEXP seqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type seq(seqSEXP);
    rcpp_result_gen = Rcpp::wrap(getFragmentIonIndices(seq));
    return rcpp_result_gen;
END_RCPP
}
// getFragmentIonSequences
Rcpp::CharacterVector getFragmentIonSequences(std::string seq);
RcppExport SEXP _RPeptideUtils_getFragmentIonSequences(SEXP seqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type seq(seqSEXP);
    rcpp_result_gen = Rcpp::wrap(getFragmentIonSequences(seq));
    return rcpp_result_gen;
END_RCPP
}
// getSequences
Rcpp::CharacterVector getSequences(const Rcpp::CharacterVector& ids, std::string fastaPath);
RcppExport SEXP _RPeptideUtils_getSequences(SEXP idsSEXP, SEXP fastaPathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ids(idsSEXP);
    Rcpp::traits::input_parameter< std::string >::type fastaPath(fastaPathSEXP);
    rcpp_result_gen = Rcpp::wrap(getSequences(ids, fastaPath));
    return rcpp_result_gen;
END_RCPP
}
// nBefore
Rcpp::CharacterVector nBefore(const Rcpp::CharacterVector& query, const Rcpp::CharacterVector ref, const Rcpp::IntegerVector& n, bool noExcept);
RcppExport SEXP _RPeptideUtils_nBefore(SEXP querySEXP, SEXP refSEXP, SEXP nSEXP, SEXP noExceptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type query(querySEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector >::type ref(refSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type n(nSEXP);
    Rcpp::traits::input_parameter< bool >::type noExcept(noExceptSEXP);
    rcpp_result_gen = Rcpp::wrap(nBefore(query, ref, n, noExcept));
    return rcpp_result_gen;
END_RCPP
}
// nAfter
Rcpp::CharacterVector nAfter(const Rcpp::CharacterVector& query, const Rcpp::CharacterVector ref, const Rcpp::IntegerVector& n, bool noExcept);
RcppExport SEXP _RPeptideUtils_nAfter(SEXP querySEXP, SEXP refSEXP, SEXP nSEXP, SEXP noExceptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type query(querySEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector >::type ref(refSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type n(nSEXP);
    Rcpp::traits::input_parameter< bool >::type noExcept(noExceptSEXP);
    rcpp_result_gen = Rcpp::wrap(nAfter(query, ref, n, noExcept));
    return rcpp_result_gen;
END_RCPP
}
// indexN
Rcpp::IntegerVector indexN(const Rcpp::CharacterVector& query, const Rcpp::CharacterVector ref, long n, bool noExcept);
RcppExport SEXP _RPeptideUtils_indexN(SEXP querySEXP, SEXP refSEXP, SEXP nSEXP, SEXP noExceptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type query(querySEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector >::type ref(refSEXP);
    Rcpp::traits::input_parameter< long >::type n(nSEXP);
    Rcpp::traits::input_parameter< bool >::type noExcept(noExceptSEXP);
    rcpp_result_gen = Rcpp::wrap(indexN(query, ref, n, noExcept));
    return rcpp_result_gen;
END_RCPP
}
// getModifiedResidues
Rcpp::CharacterVector getModifiedResidues(const Rcpp::CharacterVector& ids, const Rcpp::CharacterVector& peptideSeq, std::string fastaPath, std::string modSep);
RcppExport SEXP _RPeptideUtils_getModifiedResidues(SEXP idsSEXP, SEXP peptideSeqSEXP, SEXP fastaPathSEXP, SEXP modSepSEXP) {
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
RcppExport SEXP _RPeptideUtils_combineMods(SEXP modsSEXP, SEXP sepSEXP) {
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
Rcpp::NumericVector calcMass(const Rcpp::StringVector& sequences, bool monoMass, std::string residueAtoms);
RcppExport SEXP _RPeptideUtils_calcMass(SEXP sequencesSEXP, SEXP monoMassSEXP, SEXP residueAtomsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type sequences(sequencesSEXP);
    Rcpp::traits::input_parameter< bool >::type monoMass(monoMassSEXP);
    Rcpp::traits::input_parameter< std::string >::type residueAtoms(residueAtomsSEXP);
    rcpp_result_gen = Rcpp::wrap(calcMass(sequences, monoMass, residueAtoms));
    return rcpp_result_gen;
END_RCPP
}
// calcFormula
Rcpp::StringVector calcFormula(const Rcpp::StringVector& sequences, bool subscripts, std::string residueAtoms);
RcppExport SEXP _RPeptideUtils_calcFormula(SEXP sequencesSEXP, SEXP subscriptsSEXP, SEXP residueAtomsSEXP) {
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
RcppExport SEXP _RPeptideUtils_oneLetterToThree(SEXP sequencesSEXP, SEXP sep_inSEXP, SEXP sep_outSEXP, SEXP n_term_outSEXP, SEXP c_term_outSEXP) {
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
RcppExport SEXP _RPeptideUtils_threeLetterToOne(SEXP sequencesSEXP, SEXP sep_inSEXP, SEXP sep_outSEXP, SEXP n_term_outSEXP, SEXP c_term_outSEXP) {
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
// readFasta
Rcpp::DataFrame readFasta(std::string fastaPath, long n_entries);
RcppExport SEXP _RPeptideUtils_readFasta(SEXP fastaPathSEXP, SEXP n_entriesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fastaPath(fastaPathSEXP);
    Rcpp::traits::input_parameter< long >::type n_entries(n_entriesSEXP);
    rcpp_result_gen = Rcpp::wrap(readFasta(fastaPath, n_entries));
    return rcpp_result_gen;
END_RCPP
}
// fastaInfo
Rcpp::List fastaInfo(std::string fastaPath);
RcppExport SEXP _RPeptideUtils_fastaInfo(SEXP fastaPathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fastaPath(fastaPathSEXP);
    rcpp_result_gen = Rcpp::wrap(fastaInfo(fastaPath));
    return rcpp_result_gen;
END_RCPP
}
// transpose_sequence
Rcpp::DataFrame transpose_sequence(const Rcpp::StringVector& peptide_sequences, const Rcpp::NumericVector& quantification, const std::string& protein_seq);
RcppExport SEXP _RPeptideUtils_transpose_sequence(SEXP peptide_sequencesSEXP, SEXP quantificationSEXP, SEXP protein_seqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type peptide_sequences(peptide_sequencesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type quantification(quantificationSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type protein_seq(protein_seqSEXP);
    rcpp_result_gen = Rcpp::wrap(transpose_sequence(peptide_sequences, quantification, protein_seq));
    return rcpp_result_gen;
END_RCPP
}
// digest
Rcpp::List digest(Rcpp::CharacterVector sequences, Rcpp::CharacterVector ids, unsigned nMissedCleavages, std::string cleavagePattern, bool mz_filter, std::string residueAtoms, double minMz, double maxMz, int minCharge, int maxCharge, size_t minLen, size_t maxLen);
RcppExport SEXP _RPeptideUtils_digest(SEXP sequencesSEXP, SEXP idsSEXP, SEXP nMissedCleavagesSEXP, SEXP cleavagePatternSEXP, SEXP mz_filterSEXP, SEXP residueAtomsSEXP, SEXP minMzSEXP, SEXP maxMzSEXP, SEXP minChargeSEXP, SEXP maxChargeSEXP, SEXP minLenSEXP, SEXP maxLenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type sequences(sequencesSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type ids(idsSEXP);
    Rcpp::traits::input_parameter< unsigned >::type nMissedCleavages(nMissedCleavagesSEXP);
    Rcpp::traits::input_parameter< std::string >::type cleavagePattern(cleavagePatternSEXP);
    Rcpp::traits::input_parameter< bool >::type mz_filter(mz_filterSEXP);
    Rcpp::traits::input_parameter< std::string >::type residueAtoms(residueAtomsSEXP);
    Rcpp::traits::input_parameter< double >::type minMz(minMzSEXP);
    Rcpp::traits::input_parameter< double >::type maxMz(maxMzSEXP);
    Rcpp::traits::input_parameter< int >::type minCharge(minChargeSEXP);
    Rcpp::traits::input_parameter< int >::type maxCharge(maxChargeSEXP);
    Rcpp::traits::input_parameter< size_t >::type minLen(minLenSEXP);
    Rcpp::traits::input_parameter< size_t >::type maxLen(maxLenSEXP);
    rcpp_result_gen = Rcpp::wrap(digest(sequences, ids, nMissedCleavages, cleavagePattern, mz_filter, residueAtoms, minMz, maxMz, minCharge, maxCharge, minLen, maxLen));
    return rcpp_result_gen;
END_RCPP
}
// spaceCoveragePlot
Rcpp::IntegerVector spaceCoveragePlot(const Rcpp::IntegerVector& begin_i, const Rcpp::IntegerVector& end_i);
RcppExport SEXP _RPeptideUtils_spaceCoveragePlot(SEXP begin_iSEXP, SEXP end_iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type begin_i(begin_iSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type end_i(end_iSEXP);
    rcpp_result_gen = Rcpp::wrap(spaceCoveragePlot(begin_i, end_i));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RPeptideUtils_getFragmentIonIndices", (DL_FUNC) &_RPeptideUtils_getFragmentIonIndices, 1},
    {"_RPeptideUtils_getFragmentIonSequences", (DL_FUNC) &_RPeptideUtils_getFragmentIonSequences, 1},
    {"_RPeptideUtils_getSequences", (DL_FUNC) &_RPeptideUtils_getSequences, 2},
    {"_RPeptideUtils_nBefore", (DL_FUNC) &_RPeptideUtils_nBefore, 4},
    {"_RPeptideUtils_nAfter", (DL_FUNC) &_RPeptideUtils_nAfter, 4},
    {"_RPeptideUtils_indexN", (DL_FUNC) &_RPeptideUtils_indexN, 4},
    {"_RPeptideUtils_getModifiedResidues", (DL_FUNC) &_RPeptideUtils_getModifiedResidues, 4},
    {"_RPeptideUtils_combineMods", (DL_FUNC) &_RPeptideUtils_combineMods, 2},
    {"_RPeptideUtils_calcMass", (DL_FUNC) &_RPeptideUtils_calcMass, 3},
    {"_RPeptideUtils_calcFormula", (DL_FUNC) &_RPeptideUtils_calcFormula, 3},
    {"_RPeptideUtils_oneLetterToThree", (DL_FUNC) &_RPeptideUtils_oneLetterToThree, 5},
    {"_RPeptideUtils_threeLetterToOne", (DL_FUNC) &_RPeptideUtils_threeLetterToOne, 5},
    {"_RPeptideUtils_readFasta", (DL_FUNC) &_RPeptideUtils_readFasta, 2},
    {"_RPeptideUtils_fastaInfo", (DL_FUNC) &_RPeptideUtils_fastaInfo, 1},
    {"_RPeptideUtils_transpose_sequence", (DL_FUNC) &_RPeptideUtils_transpose_sequence, 3},
    {"_RPeptideUtils_digest", (DL_FUNC) &_RPeptideUtils_digest, 12},
    {"_RPeptideUtils_spaceCoveragePlot", (DL_FUNC) &_RPeptideUtils_spaceCoveragePlot, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_RPeptideUtils(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

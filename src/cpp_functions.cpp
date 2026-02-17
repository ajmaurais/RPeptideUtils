
#include <Rcpp.h>
#include <stdexcept>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <thread>
#include <atomic>
#include <chrono>
#include <regex>

#include <cmath>
#include <fstream>
#include <sstream>

#include <fastaFile.hpp>
#include <molecularFormula.hpp>
#include <sequenceUtils.hpp>
#include <utils.hpp>
#include <msInterface/mzMLFile.hpp>

//!Return data files included in an R package
std::string _getPackageData(std::string filename, std::string packageName = "RPeptideUtils")
{
    Rcpp::Environment base("package:base");
    Rcpp::Function sys_file = base["system.file"];
    Rcpp::StringVector file_path_sv = sys_file(
        filename,
        Rcpp::_["package"] = packageName,
        Rcpp::_["mustWork"] = true
    );
    std::string file_path = Rcpp::as<std::string>(file_path_sv);
    return file_path;
}

//' Get indices of fragment ions for peptide sequence.
//'
//' @title Get fragment ion sequence indices.
//' @param seq Peptide sequence as a string.
//' @return fragments IntegerMatrix where each row is a fragment ion. The first column is the startig index of the fragment ion in the peptide sequence and the second column is the fragment length.
// [[Rcpp::export]]
Rcpp::IntegerMatrix getFragmentIonIndices(std::string seq)
{
    std::map<std::string, utils::SizePair> ions;
    utils::seqToIons(seq, ions);
    Rcpp::IntegerMatrix fragments(ions.size(), 2);
    Rcpp::CharacterVector keys;
    Rcpp::IntegerVector start, length;
    for(auto ion: ions){
        keys.push_back(ion.first.c_str());
        start.push_back(ion.second.first);
        length.push_back(ion.second.second);
    }
    rownames(fragments) = keys;
    colnames(fragments) = Rcpp::CharacterVector::create("start", "length");
    fragments(Rcpp::_, 0) = start;
    fragments(Rcpp::_, 1) = length;
    return fragments;
}

//' Get sequences of fragment ions for peptide sequence.
//'
//' @title Get fragment ion sequences.
//' @param seq Peptide sequence as a string.
//' @return fragments Named CharacterVector where the names are the ions and the values are the sequences.
// [[Rcpp::export]]
Rcpp::CharacterVector getFragmentIonSequences(std::string seq)
{
    std::map<std::string, utils::SizePair> ions;
    utils::seqToIons(seq, ions);
    Rcpp::CharacterVector keys;
    Rcpp::CharacterVector fragments;
    for(auto ion: ions){
        keys.push_back(ion.first.c_str());
        fragments.push_back(seq.substr(ion.second.first, ion.second.second).c_str());
    }
    fragments.names() = keys;
    return fragments;
}

//' Get protein sequences for a vector of uniprot IDs
//'
//' @title Get protein IDs in fasta file for a vector of Uniprot IDs
//' @param ids CharacterVector of uniprot IDs
//' @param fastaPath path to fasta formatted file to look up protein sequences.
//' @return CharacterVector of protein sequences in same order as ids.
//'
//' @examples
//' #By default the fasta file included in the package containing human protein sequences is used.
//' getSequences(c("A0MZ66", "A6NMY6", "O00213", "O00213"))
//'
//' #A fasta file can also be manually specified.
//' fasta_path <- system.file('extdata/Human_uniprot-reviewed_20171020.fasta', package = 'RPeptideUtils')
//' getSequences(c("A0MZ66", "A6NMY6", "O00213", "O00213"), fasta_path)
//'
// [[Rcpp::export]]
Rcpp::CharacterVector getSequences(const Rcpp::CharacterVector& ids, std::string fastaPath = "")
{
    std::string _fastaPath = fastaPath.empty() ? _getPackageData("extdata/Human_uniprot-reviewed_20171020.fasta") : fastaPath;

    Rcpp::CharacterVector ret;

    utils::FastaFile fasta(true, utils::absPath(_fastaPath));
    if(!fasta.read()) throw std::runtime_error("Could not read fasta file!");

    size_t len = ids.size();
    for(size_t i = 0; i < len; i++){
        ret.push_back(fasta.getSequence(std::string(ids[i])));
    }

    return ret;
}

//' Get n residues before query in ref. If n overruns ref, the maximum number of characters will be returned.
//'
//' @title get n residues before query.
//' @param query String to search for.
//' @param ref String to search in.
//' @param n Number of residues in output.
//' @param noExcept Should an std::runtime_error be thrown if query is not in ref?
//'
//' @return n residues before query.
//'
// [[Rcpp::export]]
Rcpp::CharacterVector nBefore(const Rcpp::CharacterVector& query, const Rcpp::CharacterVector ref,
                              size_t n, bool noExcept = false)
{
    Rcpp::CharacterVector ret;

    size_t len = ref.size();
    if(len != 1){
        if(len != query.size())
            throw std::runtime_error("query and ref must be the same length!");
    }

    for(size_t i = 0; i < query.size(); i++){
        std::string ref_temp = std::string(len == 1 ? ref[0] : ref[i]);
        ret.push_back(utils::nBefore(std::string(query[i]), ref_temp, n, noExcept));
    }

    return ret;
}

//' Get n residues after query in ref. If n overruns ref, the maximum number of characters will be returned.
//'
//' @title get n residues after query.
//' @param query String to search for.
//' @param ref String to search in.
//' @param n Number of residues in output.
//' @param noExcept Should an std::runtime_error be thrown if query is not in ref?
//'
//' @return n residues after query.
//'
// [[Rcpp::export]]
Rcpp::CharacterVector nAfter(const Rcpp::CharacterVector& query, const Rcpp::CharacterVector ref,
                             size_t n, bool noExcept = false)
{
    Rcpp::CharacterVector ret;

    size_t ref_len = ref.size();
    if(ref_len != 1){
        if(ref_len != query.size())
            throw std::runtime_error("query and ref must be the same length!");
    }

    for(size_t i = 0; i < query.size(); i++){
        std::string ref_temp = std::string(ref_len == 1 ? ref[0] : ref[i]);
        ret.push_back(utils::nAfter(std::string(query[i]), ref_temp, n, noExcept));
    }

    return ret;
}

//' Get the index of residue n of query in ref. If n is -1, the index of the last residue in query is returned.
//'
//' @title get the index of nth residues of query in ref.
//' @param query String to search for.
//' @param ref String to search in.
//' @param n Residue number in query.
//' @param noExcept Should an std::runtime_error be thrown if query is not in ref?
//'
//' @return n residues after query.
//'
// [[Rcpp::export]]
Rcpp::IntegerVector indexN(const Rcpp::CharacterVector& query, const Rcpp::CharacterVector ref,
    long n = 1, bool noExcept = false)
{
    Rcpp::IntegerVector ret;

    size_t len = ref.size();
    if(len != 1){
        if(len != query.size())
            throw std::runtime_error("query and ref must be the same length!");
    }

    for(size_t i = 0; i < query.size(); i++){
        std::string ref_temp = std::string(len == 1 ? ref[0] : ref[i]);
        ret.push_back(utils::indexN(std::string(query[i]), ref_temp, n == -1 ? std::string::npos : n, noExcept));
    }

    return ret;
}

//remove residues before and after cleavage
std::string makeSequenceFromFullSequence(std::string fs)
{
    fs = fs.substr(fs.find(".") + 1);
    fs = fs.substr(0, fs.find_last_of("."));
    return fs;
}

//' Get locations of modified residues in parent protein
//'
//' @title Get locations of modified residues in parent protein
//' @param fastaPath path to fasta formatted file to look up protein sequences
//' @param ids CharacterVector of Uniprot IDs
//' @param peptideSeq CharacterVector of peptide sequences containing modifications
//' @param modSep delimiter for multiple modifications
//' @return CharacterVector containing locations of modifications in protein sequence
//'
//' @examples
//' getModifiedResidues(c("Q00839", "Q9HCS7", "Q7L014"), c("APQC*LGK", "FADMEC*K", "GAEIIVC*TPGR"))
//'
// [[Rcpp::export]]
Rcpp::CharacterVector getModifiedResidues(const Rcpp::CharacterVector& ids,
    const Rcpp::CharacterVector& peptideSeq, std::string fastaPath = "", std::string modSep = "|")
{
    std::string _fastaPath = fastaPath.empty() ? _getPackageData("extdata/Human_uniprot-reviewed_20171020.fasta") : fastaPath;

    size_t len = ids.size();
    if(len != peptideSeq.size())
        throw std::runtime_error("ids.size() != peptideSeq.size()");

    //init FastaFile
    utils::FastaFile fasta(true, _fastaPath);
        if(!fasta.read()) throw std::runtime_error("Could not read fasta file!");

    Rcpp::CharacterVector ret(len, "");
    std::vector<int> modIndex_temp;
    for(size_t i = 0; i < len; i++)
    {
        std::string this_modLocs;
        std::string seqTemp = utils::getModLocs(std::string(peptideSeq[i]), modIndex_temp);
        for(auto it = modIndex_temp.begin(); it != modIndex_temp.end(); ++it){
            utils::addChar(fasta.getModifiedResidue(std::string(ids[i]), seqTemp, *it), this_modLocs, modSep);
        }
        ret[i] = this_modLocs;
    }
    return ret;
}

//' Combine concated mods from multiple peptides into a single string.
//'
//' @title Combine mods from multiple peptides into a single string
//' @param mods Modifications to combine
//' @param sep delimiter separating modifications
//' @return Modifications combined into a single string
//'
//' @examples
//' combineMods(c('C157', 'C157|C125', 'C50', 'C125'))
//'
// [[Rcpp::export]]
std::string combineMods(const Rcpp::CharacterVector& mods, char sep = '|')
{
    std::set<std::string> found;

    size_t len = mods.size();
    for(size_t i = 0; i < len; i++)
    {
        std::vector<std::string> temp;
        utils::split(std::string(mods[i]), sep, temp);
        found.insert(temp.begin(), temp.end());
    }

    return utils::concat(found.begin(), found.end());
}

//' Calculate peptide monoisotopic or average masses.
//'
//' @title Calculate peptide masses
//' @param sequences Peptide sequences
//' @param monoMass Should monoisotopic mass be calculated. If false, average mass is calculated.
//' @param residueAtoms Path to residueAtoms file. If blank, the default file included in the package is used.
//' @return vector of peptide masses.
//'
//' @examples
//' calcMass(c("ACLLPETVNMEEYPYDAEY", "ALCAEFK"))
//'
// [[Rcpp::export]]
Rcpp::NumericVector calcMass(const Rcpp::StringVector& sequences,
                             bool monoMass = true,
                             std::string residueAtoms = "")
{
    //get data file paths
    std::string residueAtomsPath = residueAtoms.empty() ? _getPackageData("defaultResidueAtoms.txt") : residueAtoms;
    char avg_mono = monoMass ? 'm' : 'a';

    //init residues
    utils::Residues residues(residueAtomsPath);
    if(!residues.initialize(false)) throw std::runtime_error("Error reading required files for calcMass!");

    size_t len = sequences.size();
    Rcpp::NumericVector ret(len);
    for(size_t i = 0; i < len; i++){
        ret[i] = residues.calcMass(std::string(sequences[i]), avg_mono);
    }

    return ret;
}

//' Calculate peptide molecular formulas
//'
//' @title Calculate peptide formulas
//' @param sequences Peptide sequences
//' @param subscripts Should formulas have subscripts or normal baseline numbers?
//' @param residueAtoms Path to residueAtoms file. If blank, the default file included in the package is used.
//' @return vector of peptide formulas.
//'
//' @examples
//' calcFormula(c("ACLLPETVNMEEYPYDAEY", "ALCAEFK"), subscripts = TRUE)
//'
// [[Rcpp::export]]
Rcpp::StringVector calcFormula(const Rcpp::StringVector& sequences,
                                                             bool subscripts = false,
                                                             std::string residueAtoms = "")
{
    //get data file paths
    std::string residueAtomsPath = residueAtoms.empty() ? _getPackageData("defaultResidueAtoms.txt") : residueAtoms;

    //init residues
    utils::Residues residues;
    if(!residues.initialize(residueAtomsPath))
        throw std::runtime_error("Error reading required files for calcFormula!");

    size_t len = sequences.size();
    Rcpp::StringVector ret(len);
    for(size_t i = 0; i < len; i++){
        ret[i] = residues.calcFormula(std::string(sequences[i]), subscripts);
    }

    return ret;
}

//' Convert from 1 letter amino acid codes to 3
//'
//' @title Convert to 3 letter amino acid codes
//' @param sequences vector of sequences
//' @param sep_in deliminator between amino acids in input
//' @param sep_out deliminator between amino acids in output
//' @param n_term_out string to append to n terminus
//' @param c_term_out string to append to c terminus
//' @return StringVector of peptides with three letter amino acid codes
//'
//' @examples
//' oneLetterToThree(c("AC*LLPETVNMEEYPYDAEY", "ALCAEFK", "AQUPIVER", "C*TGGEVGATSALAPK"))
//'
// [[Rcpp::export]]
Rcpp::StringVector oneLetterToThree(Rcpp::StringVector sequences,
                                    std::string sep_in = "",
                                    std::string sep_out = "",
                                    std::string n_term_out = "",
                                    std::string c_term_out = "")
{
    size_t len = sequences.size();
    Rcpp::StringVector ret(len);
    for(size_t i = 0; i < len; i++){
        ret[i] = utils::oneLetterToThree(std::string(sequences[i]),
                                                                         sep_in, sep_out,
                                                                         n_term_out, c_term_out);
    }
    return ret;
}

//' Convert from 3 letter amino acid codes to 1
//'
//' @title Convert to 1 letter amino acid codes
//' @param sequences vector of sequences
//' @param sep_in deliminator between amino acids in input
//' @param sep_out deliminator between amino acids in output
//' @param n_term_out string to append to n terminus
//' @param c_term_out string to append to c terminus
//' @return StringVector of peptides with one letter amino acid codes
//'
//' @examples
//' threeLetterToOne(c("Ala-Cys*-Leu-Leu-Pro", "Ala-Leu-Cys-Ala", "Ala-Gln-Sec-Ile"), sep_in = "-")
//'
// [[Rcpp::export]]
Rcpp::StringVector threeLetterToOne(Rcpp::StringVector sequences,
                                    std::string sep_in = "",
                                    std::string sep_out = "",
                                    std::string n_term_out = "",
                                    std::string c_term_out = "")
{
    size_t len = sequences.size();
    Rcpp::StringVector ret(len);
    for(size_t i = 0; i < len; i++){
        ret[i] = utils::threeLetterToOne(std::string(sequences[i]),
                                                                         sep_in, sep_out,
                                                                         n_term_out, c_term_out);
    }
    return ret;
}


//' Read all sequences in fasta file. Reverse matches are automatically skipped.
//'
//' @title Read fasta file.
//'
//' @param fastaPath Path to fasta file. Be default, fasta file included in package is used.
//' @param n_entries Number of entries to read. If 0, all entries are read.
//' @return DataFrame with columns for ID and sequence.
//'
// [[Rcpp::export]]
Rcpp::DataFrame readFasta(std::string fastaPath = "", long n_entries = 0)
{
    std::string _fastaPath = fastaPath.empty() ?
        _getPackageData("extdata/Human_uniprot-reviewed_20171020.fasta") : fastaPath;

    //init FastaFile
    utils::FastaFile fasta(true, _fastaPath);
        if(!fasta.read()) throw std::runtime_error("Could not read fasta file!");

    size_t totalEntries = fasta.getSequenceCount();
    if(n_entries < 0 || static_cast<size_t>(n_entries) > totalEntries)
        throw std::runtime_error("n_entries more than the number of entries in file!");
    size_t len = n_entries == 0 ? totalEntries : static_cast<size_t>(n_entries);

    // Pre-allocate vectors to avoid repeated reallocations
    Rcpp::CharacterVector ids(len);
    Rcpp::CharacterVector seqs(len);

    for(size_t i = 0; i < len; i++)
    {
        ids[i] = fasta.getIndexID(i);
        seqs[i] = fasta.at(i);
    }

    return Rcpp::DataFrame::create(Rcpp::Named("id") = ids,
                                   Rcpp::Named("sequence") = seqs,
                                   Rcpp::Named("stringsAsFactors") = false);
}

//' Get metadata about a fasta file.
//'
//' @title Get fasta file info.
//'
//' @param fastaPath Path to fasta file. Be default, fasta file included in package is used.
//' @return List with slots for sequence count and vector of entry IDs contained in file.
//'
// [[Rcpp::export]]
Rcpp::List fastaInfo(std::string fastaPath = "")
{
    std::string _fastaPath = fastaPath.empty() ?
        _getPackageData("extdata/Human_uniprot-reviewed_20171020.fasta") : fastaPath;

    //init FastaFile
    utils::FastaFile fasta(true, _fastaPath);
        if(!fasta.read()) throw std::runtime_error("Could not read fasta file!");

    size_t len = fasta.getSequenceCount();
    Rcpp::CharacterVector ids(len, "");
    for(size_t i = 0; i < len; i++){
        ids[i] = fasta.getIndexID(i);
    }

    return Rcpp::List::create(Rcpp::Named("seq_count") = len,
                              Rcpp::Named("path") = fastaPath,
                              Rcpp::Named("ids") = ids);
}

//' Transpose peptide quantifications for a single protein into amino acid level
//' quantifications. A row will be included in the output for each time an amino acid
//' at a given position was included in a peptide in peptide_sequences. Additional processing
//' is required to obtain summary values for each amino acid position.
//'
//' @title Transpose peptide quantifications into amino acid level quantifications.
//'
//' @param peptide_sequences List of peptide sequences.
//' @param quantification Ratio or spectral count values for peptide_sequences.
//' @param protein_seq Parent protein sequence.
//' @return DataFrame with columns for 'residue', 'number', and 'quant'.
//'
// [[Rcpp::export]]
Rcpp::DataFrame transpose_sequence(const Rcpp::StringVector& peptide_sequences,
                                   const Rcpp::NumericVector& quantification,
                                   const std::string& protein_seq)
{
    if(peptide_sequences.size() != quantification.size())
        throw std::runtime_error("peptide_sequences and quantification must be the same length!");

    std::vector<char> residues;
    Rcpp::IntegerVector numbers;
    Rcpp::NumericVector quantifications;

    size_t n_seq = peptide_sequences.size();
    size_t begin, end;
    for(size_t i = 0; i < n_seq; i++)
    {
        if(!utils::align(std::string(peptide_sequences[i]), protein_seq, begin, end))
            throw std::runtime_error("Peptide sequence '" + peptide_sequences[i] + "' does not exist in protein_seq!");

        size_t pep_len = peptide_sequences[i].size();
        for(size_t pep_begin = 0; pep_begin < pep_len; pep_begin++) {
            residues.push_back(peptide_sequences[i][pep_begin]);
            numbers.push_back(pep_begin + begin);
            quantifications.push_back(quantification[i]);
        }
    }

    return Rcpp::DataFrame::create(Rcpp::Named("residue") = residues,
                                   Rcpp::Named("number") = numbers,
                                   Rcpp::Named("quant") = quantifications,
                                   Rcpp::Named("stringsAsFactors") = false);
}


//' The function uses charge and m/z filters to remove peptides which would not be
//' observable by MS. The m/z for peptides in charge states minCharge to maxCharge
//' are calculated. If the m/z for any charge state is in between minMZ and maxMZ, the
//' sequence will be appended to peptides.
//'
//' @title Perform a virtual protease digest of a protein.
//'
//' @param sequences StringVector containing protein sequences. Whitespace will automatically be removed.
//' @param ids Names for the slot for each protein's peptides in output.
//' @param nMissedCleavages number of missed cleavages to allow.
//' @param cleavagePattern RegEx for protease cleavage pattern. The default is the pattern for trypsin.
//' @param mz_filter Should peptides included in output be filtered by mz?
//' @param residueAtoms Path to residueAtoms file. If blank, the default file included in the package is used.
//' @param minMz Minimum m/z to allow in peptides.
//' @param maxMz Maximum m/z to allow in peptides. Set to 0 for no upper bound on m/z.
//' @param minCharge Minimum charge to consider when calculating m/z.
//' @param maxCharge Maximum charge to consider when calculating m/z.
//' @param minLen Minimum peptide length.
//' @param maxLen Maximum peptide length. Set to 0 for no upper bound on length.
//'
//' @return A list with named elements containing vectors of each input protein's peptides.
//'
//' @examples
//' digest(c("KLGAARKLGAGLAKVIGAGIGIGK", "KLGAARKLGAGLAKPVIGAGIGIGK"), c('a', 'b'))
//'
// [[Rcpp::export]]
Rcpp::List digest(Rcpp::CharacterVector sequences, Rcpp::CharacterVector ids,
                  unsigned nMissedCleavages = 0, std::string cleavagePattern = "([RK])(?=[^P])",
                  bool mz_filter = true, std::string residueAtoms = "",
                  double minMz = 400, double maxMz = 1800,
                  int minCharge = 1, int maxCharge = 5,
                  size_t minLen = 6, size_t maxLen = 0)
{
    //check args
    if(ids.size() != sequences.size())
        throw std::runtime_error("sequences and ids must be the same length!");

    size_t _maxLen = maxLen == 0 ? std::string::npos : maxLen;
    Rcpp::List ret;
    size_t len = sequences.size();

    // Compile regex once outside the loop
    std::regex re(cleavagePattern);

    // Initialize Residues once outside the loop (only if mz_filter is true)
    utils::Residues residues;
    if(mz_filter) {
        std::string residueAtomsPath = residueAtoms.empty() ?
            _getPackageData("defaultResidueAtoms.txt") : residueAtoms;
        residues = utils::Residues(residueAtomsPath);
        if(!residues.initialize(false))
            throw std::runtime_error("Error reading required files for calcMass!");
    }

    for(size_t i = 0; i < len; i++)
    {
        std::vector<std::string> peptides_temp;
        Rcpp::CharacterVector ret_temp;
        std::string seq = utils::removeWhitespace(std::string(sequences[i]));

        if(mz_filter)
        {
            residues.digest(seq, peptides_temp, nMissedCleavages, false, re,
                            minMz, maxMz, minCharge, maxCharge);

            std::sort(peptides_temp.begin(), peptides_temp.end(), utils::strLenCompare());
            for(const auto& pep : peptides_temp)
            {
                size_t pep_len = pep.length();
                if(pep_len >= minLen && (_maxLen == std::string::npos || pep_len <= _maxLen))
                    ret_temp.push_back(pep.c_str());
            }
        }
        else {
            utils::digest(seq, peptides_temp, nMissedCleavages, minLen, _maxLen, re);
            for(const auto& pep : peptides_temp) {
                ret_temp.push_back(pep.c_str());
            }
        }
        ret.push_back(ret_temp, std::string(ids[i]));
    }

    return ret;
}

void matchingProteinsWorker(const std::map<std::string, std::string>& proteins,
                       std::map<std::string, std::vector<std::string>  >& matchesIds,
                       std::atomic<size_t>& peptideIndex) {
    for(auto peptide: matchesIds) {
        for(auto prot: proteins) {
            if(prot.second.find(peptide.first) != std::string::npos) {
                matchesIds[peptide.first].push_back(prot.first);
            }
        }
        peptideIndex++;
    }
}

void progressBarWorker(std::atomic<size_t>& index, size_t count,
                       const std::string& message, int sleepTime)
{
    size_t lastIndex = index.load();
    size_t curIndex = lastIndex;
    int noChangeIterations = 0;

    Rcpp::Rcout << message << '\n';
    while(index < count) {
        curIndex = index.load();

        if(lastIndex == curIndex) noChangeIterations++;
        else noChangeIterations = 0;

        if(noChangeIterations > 30)
            throw std::runtime_error("Thread timeout!");

        utils::printProgress(float(curIndex) / float(count));
        std::this_thread::sleep_for(std::chrono::seconds(sleepTime));

        lastIndex = curIndex;
    }
    utils::printProgress(float(index.load()) / float(count));
    Rcpp::Rcout << NEW_LINE;
    Rcpp::Rcout << "Done!\n";
}

//' Given a vector of peptide sequences, find all the proteins containing the peptide in a fasta file.
//'
//' @title Find all proteins containing peptide sequences.
//' @param peptides A character vector of peptide sequences.
//' @param fastaPath path to fasta formatted file to look up protein sequences.
//' @param progressBar Show progress bar?
//' @param nThread Number of threads to use. By default use 1 thread per virtual core on machine.
//' @return List where names are peptide sequences, and values are the ids for the matching proteins.
//'
//' @examples
//' fasta_path <- system.file('extdata/Human_uniprot-reviewed_20171020.fasta', package = 'RPeptideUtils')
//' matchingProteins(c('PEPTIDE'), progressBar = FALSE, fastaPath = fasta_path)
//'
// [[Rcpp::export]]
Rcpp::List matchingProteins(Rcpp::CharacterVector peptides, std::string fastaPath = "",
                       bool progressBar = true, size_t nThread = 0)
{
    // convert peptide seq char*(s) to std::string
    size_t len = peptides.size();
    std::vector<std::string> peptides_s;
    for(size_t i = 0; i < len; i++) {
        peptides_s.push_back(std::string(peptides[i]));
    }

    // read fasta file
    std::string _fastaPath = fastaPath.empty() ?
        _getPackageData("extdata/Human_uniprot-reviewed_20171020.fasta") : fastaPath;

    utils::FastaFile fasta(true, _fastaPath);
    if(!fasta.read()) throw std::runtime_error("Could not read fasta file!");
    if(progressBar) Rcpp::Rcout << "Reading fasta file..." << std::flush;
    std::map<std::string, std::string> sequences;
    size_t n_seq = fasta.getSequenceCount();
    for(size_t i = 0; i < n_seq; i++) {
        sequences[fasta.getIndexID(i)] = fasta.at(i);
    }
    if(progressBar) Rcpp::Rcout << " Done!\n";

    // split up input to fit worker threads
    if(nThread == 0)
        nThread = std::min(size_t(std::thread::hardware_concurrency()), len);
    size_t peptides_per_thread = len / nThread;
    if(len % nThread != 0) peptides_per_thread += 1;

    // init threads
    std::vector<std::thread> threads;
    std::atomic<size_t> peptideIndex(0);
    if(progressBar) {
        std::string message = "Searching matching protein sequences for " + std::to_string(len) +
                              " peptides using " + std::to_string(nThread) + " threads...";
        threads.emplace_back(progressBarWorker, std::ref(peptideIndex),
                             len, message, 1);
    }

    auto* splitPeptides = new std::map<std::string, std::vector<std::string> >[nThread];
    size_t begin, end;
    for(size_t i = 0; i < nThread; i++) {
        begin = i * peptides_per_thread;
        end = (begin + peptides_per_thread > len ? len : begin + peptides_per_thread);

        for(size_t j = begin; j < end; j++){
            splitPeptides[i][peptides_s.at(j)] = std::vector<std::string> ();
        }
        threads.emplace_back(matchingProteinsWorker, std::ref(sequences),
                             std::ref(splitPeptides[i]),
                             std::ref(peptideIndex));
    }

    // join threads
    for(auto& t : threads) {
        t.join();
    }

    if(progressBar) Rcpp::Rcout << "Generating output...\n";
    Rcpp::List ret;
    for(size_t i = 0; i < nThread; i++) {
        for(auto& it : splitPeptides[i]) {
            if(progressBar) {
                if(ret.size() % 1000 == 0) {
                    utils::printProgress(float(ret.size()) / float(len));
                }
            }

            ret.push_back(it.second, it.first);
        }
    }
    if(progressBar){
        utils::printProgress(float(ret.size()) / float(len));
        Rcpp::Rcout << "\nDone!\n";
    }

    delete [] splitPeptides;
    return ret;
}

struct ScanResult {
    size_t originalIndex;
    size_t scanNum;
    std::string fileBasename;
    std::string filePath;
    double rt;
    int msLevel;
    double isoWinTarget;
    double isoWinLower;
    double isoWinUpper;
};

struct FileGroup {
    std::string filePath;
    std::vector<std::pair<size_t, size_t>> tasks; // (originalIndex, scanNum)
};

//! Extract a double value attribute from a cvParam XML string
double _extractCvParamValue(const std::string& xml, const std::string& accession)
{
    size_t pos = xml.find(accession);
    if(pos == std::string::npos) return std::numeric_limits<double>::quiet_NaN();
    size_t valPos = xml.find("value=\"", pos);
    if(valPos == std::string::npos) return std::numeric_limits<double>::quiet_NaN();
    valPos += 7;
    size_t valEnd = xml.find("\"", valPos);
    if(valEnd == std::string::npos) return std::numeric_limits<double>::quiet_NaN();
    return std::stod(xml.substr(valPos, valEnd - valPos));
}

//! Parse isolation windows from raw mzML file content.
//! Returns a map from scan number to (targetMz, lowerOffset, upperOffset).
std::map<size_t, std::tuple<double, double, double>>
_parseIsolationWindows(const std::string& content)
{
    std::map<size_t, std::tuple<double, double, double>> result;
    size_t pos = 0;

    while((pos = content.find("<spectrum ", pos)) != std::string::npos) {
        size_t endPos = content.find("</spectrum>", pos);
        if(endPos == std::string::npos) break;

        // Extract scan number from id="... scan=N ..."
        size_t idPos = content.find("id=\"", pos);
        if(idPos == std::string::npos || idPos > endPos) { pos = endPos; continue; }
        size_t scanPos = content.find("scan=", idPos);
        if(scanPos == std::string::npos || scanPos > endPos) { pos = endPos; continue; }
        size_t numStart = scanPos + 5;
        size_t numEnd = numStart;
        while(numEnd < content.size() && std::isdigit(content[numEnd])) numEnd++;
        if(numStart == numEnd) { pos = endPos; continue; }
        size_t scanNum = std::stoul(content.substr(numStart, numEnd - numStart));

        // Find isolationWindow within this spectrum
        size_t isoPos = content.find("<isolationWindow>", pos);
        if(isoPos != std::string::npos && isoPos < endPos) {
            size_t isoEnd = content.find("</isolationWindow>", isoPos);
            if(isoEnd != std::string::npos && isoEnd < endPos) {
                std::string isoXml = content.substr(isoPos, isoEnd - isoPos);
                double target = _extractCvParamValue(isoXml, "MS:1000827");
                double lower = _extractCvParamValue(isoXml, "MS:1000828");
                double upper = _extractCvParamValue(isoXml, "MS:1000829");
                result[scanNum] = std::make_tuple(target, lower, upper);
            }
        }

        pos = endPos;
    }

    return result;
}

void getMSScanMetadataWorker(
    const std::vector<FileGroup>& fileGroups,
    std::vector<ScanResult>& results,
    std::atomic<size_t>& scanCounter)
{
    for(const auto& group : fileGroups) {
        // Use MzMLFile to read and index the file
        utils::msInterface::MzMLFile mzml(group.filePath);
        if(!mzml.read())
            throw std::runtime_error("Could not read file: " + group.filePath);

        std::string basename = mzml.getParentFileBase();

        // Parse isolation windows from raw file content
        std::ifstream ifs(group.filePath);
        std::string rawContent((std::istreambuf_iterator<char>(ifs)),
                                std::istreambuf_iterator<char>());
        auto isoWindows = _parseIsolationWindows(rawContent);

        for(const auto& task : group.tasks) {
            ScanResult result;
            result.originalIndex = task.first;
            result.scanNum = task.second;
            result.filePath = group.filePath;
            result.fileBasename = basename;
            result.rt = NA_REAL;
            result.msLevel = NA_INTEGER;
            result.isoWinTarget = NA_REAL;
            result.isoWinLower = NA_REAL;
            result.isoWinUpper = NA_REAL;

            utils::msInterface::Scan scan;
            if(mzml.getScan(task.second, scan)) {
                result.rt = scan.getPrecursor().getRT();
                result.msLevel = scan.getLevel();
            }

            auto it = isoWindows.find(task.second);
            if(it != isoWindows.end()) {
                double target = std::get<0>(it->second);
                double lower = std::get<1>(it->second);
                double upper = std::get<2>(it->second);
                result.isoWinTarget = std::isnan(target) ? NA_REAL : target;
                result.isoWinLower = std::isnan(lower) ? NA_REAL : lower;
                result.isoWinUpper = std::isnan(upper) ? NA_REAL : upper;
            }

            results.push_back(result);
            scanCounter++;
        }
    }
}

//' Read scan metadata from mzML files.
//'
//' @title Read MS scan metadata from mzML files.
//' @param scans IntegerVector of scan numbers.
//' @param files CharacterVector of file paths to mzML files. Must be the same length as scans.
//' @param progressBar Show progress bar?
//' @param nThread Number of threads to use. By default use 1 thread per virtual core on machine.
//' @return DataFrame with columns for scan, file, path, rt, ms_level,
//'   precursor_mz, isolation_window_lower_offset, and isolation_window_upper_offset.
//'
//' @examples
//' mzml_path <- system.file('extdata/20240927_STEL1_Evo3_MW_15N_MS3_dilution_12_1.mzML',
//'                           package = 'RPeptideUtils')
//' getMSScanMetadata(c(1L, 2L, 3L), rep(mzml_path, 3), progressBar = FALSE)
//'
// [[Rcpp::export]]
Rcpp::DataFrame getMSScanMetadata(Rcpp::IntegerVector scans, Rcpp::CharacterVector files,
                                   bool progressBar = true, size_t nThread = 0)
{
    size_t len = scans.size();
    if(len != static_cast<size_t>(files.size()))
        throw std::runtime_error("scans and files must be the same length!");

    if(len == 0)
        return Rcpp::DataFrame::create(
            Rcpp::Named("scan") = Rcpp::IntegerVector(),
            Rcpp::Named("file") = Rcpp::CharacterVector(),
            Rcpp::Named("path") = Rcpp::CharacterVector(),
            Rcpp::Named("rt") = Rcpp::NumericVector(),
            Rcpp::Named("ms_level") = Rcpp::IntegerVector(),
            Rcpp::Named("precursor_mz") = Rcpp::NumericVector(),
            Rcpp::Named("isolation_window_lower_offset") = Rcpp::NumericVector(),
            Rcpp::Named("isolation_window_upper_offset") = Rcpp::NumericVector(),
            Rcpp::Named("stringsAsFactors") = false);

    // Group scans by file, keeping track of original indices
    std::map<std::string, std::vector<std::pair<size_t, size_t>>> fileGroupMap;
    for(size_t i = 0; i < len; i++) {
        fileGroupMap[std::string(files[i])].push_back({i, static_cast<size_t>(scans[i])});
    }

    // Convert to vector for splitting across threads
    std::vector<FileGroup> groups;
    for(auto& fg : fileGroupMap) {
        FileGroup g;
        g.filePath = fg.first;
        g.tasks = fg.second;
        groups.push_back(g);
    }

    size_t nFiles = groups.size();
    if(nThread == 0)
        nThread = std::min(size_t(std::thread::hardware_concurrency()), nFiles);
    if(nThread > nFiles) nThread = nFiles;

    size_t filesPerThread = nFiles / nThread;
    if(nFiles % nThread != 0) filesPerThread++;

    // Init threads
    std::vector<std::thread> threads;
    std::atomic<size_t> scanCounter(0);
    if(progressBar) {
        std::string message = "Reading scan metadata for " + std::to_string(len) +
                              " scans from " + std::to_string(nFiles) +
                              " files using " + std::to_string(nThread) + " threads...";
        threads.emplace_back(progressBarWorker, std::ref(scanCounter),
                             len, message, 1);
    }

    auto* threadGroups = new std::vector<FileGroup>[nThread];
    auto* threadResults = new std::vector<ScanResult>[nThread];
    for(size_t i = 0; i < nThread; i++) {
        size_t begin = i * filesPerThread;
        size_t end = std::min(begin + filesPerThread, nFiles);
        for(size_t j = begin; j < end; j++) {
            threadGroups[i].push_back(groups[j]);
        }
        threads.emplace_back(getMSScanMetadataWorker,
                             std::ref(threadGroups[i]),
                             std::ref(threadResults[i]),
                             std::ref(scanCounter));
    }

    // Join threads
    for(auto& t : threads) {
        t.join();
    }

    // Combine results from all threads
    std::vector<ScanResult> allResults;
    for(size_t i = 0; i < nThread; i++) {
        allResults.insert(allResults.end(), threadResults[i].begin(), threadResults[i].end());
    }
    delete[] threadGroups;
    delete[] threadResults;

    // Sort by original index to maintain input order
    std::sort(allResults.begin(), allResults.end(),
        [](const ScanResult& a, const ScanResult& b) {
            return a.originalIndex < b.originalIndex;
        });

    // Build output DataFrame
    Rcpp::IntegerVector outScans(len);
    Rcpp::CharacterVector outFiles(len);
    Rcpp::CharacterVector outPaths(len);
    Rcpp::NumericVector outRTs(len);
    Rcpp::IntegerVector outLevels(len);
    Rcpp::NumericVector outTargets(len);
    Rcpp::NumericVector outLowers(len);
    Rcpp::NumericVector outUppers(len);

    for(size_t i = 0; i < len; i++) {
        outScans[i] = allResults[i].scanNum;
        outFiles[i] = allResults[i].fileBasename;
        outPaths[i] = allResults[i].filePath;
        outRTs[i] = allResults[i].rt;
        outLevels[i] = allResults[i].msLevel;
        outTargets[i] = allResults[i].isoWinTarget;
        outLowers[i] = allResults[i].isoWinLower;
        outUppers[i] = allResults[i].isoWinUpper;
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("scan") = outScans,
        Rcpp::Named("file") = outFiles,
        Rcpp::Named("path") = outPaths,
        Rcpp::Named("rt") = outRTs,
        Rcpp::Named("ms_level") = outLevels,
        Rcpp::Named("precursor_mz") = outTargets,
        Rcpp::Named("isolation_window_lower_offset") = outLowers,
        Rcpp::Named("isolation_window_upper_offset") = outUppers,
        Rcpp::Named("stringsAsFactors") = false);
}

//' @title Remove the substring which is shared by the begining and end of all strings in a CharacterVector.
//' @param cStrings The CharacterVector.
//' @param verbose Print verbose output?
//' @return CharacterVector with common begining and ending removed.
//'
//' @examples
//' smallestDifferentStrings(c('SHARED_START_first_SHARED_END',
//'                            'SHARED_START_second_SHARED_END'))
//'
// [[Rcpp::export]]
Rcpp::CharacterVector smallestDifferentStrings(Rcpp::CharacterVector cStrings,
                                                 bool verbose = true)
{
    // get min string length
    std::vector<std::string> strings;
    size_t min_len = std::string::npos;
    for(size_t i = 0; i < cStrings.size(); i++) {
        strings.push_back(std::string(cStrings[i]));
        min_len = std::min(min_len, strings.back().size());
    }

    // debugging variables
    std::string commonStart, commonEnd = "";

    // remove common start
    bool done = false;
    for(size_t i = 0; i < min_len; i++){
        char c = strings[0][i];
        for(size_t j = 1; j < strings.size(); j++){
            if(c != strings[j][i]){
                done = true;
                break;
            }
        }
        if(done) break;
        commonStart += c;
    }
    for(size_t i = 0; i < strings.size(); i++){
        strings[i].erase(0, commonStart.size());
    }

    // remove common end
    done = false;
    for(size_t i = 0; i < min_len; i++){
        char c = *(strings[0].rbegin() + i);
        for(size_t j = 1; j < strings.size(); j++){
            if(c != *(strings[j].rbegin() + i)){
                done = true;
                break;
            }
        }
        if(done) break;
        commonEnd += c;
    }
    std::reverse(commonEnd.begin(), commonEnd.end());
    for(size_t i = 0; i < strings.size(); i++){
        strings[i].erase(strings[i].size() - commonEnd.size());
    }

    if(verbose)
        Rcpp::Rcout << "commonStart = " << commonStart << '\n'
                    << "commonEnd = " << commonEnd << '\n';

    Rcpp::CharacterVector ret;
    for(auto s: strings) ret.push_back(s.c_str());

    return ret;
}

//' Count the number of missed protease cleavages in peptide sequences.
//'
//' @title Count missed protease cleavages
//' @param sequences CharacterVector of peptide sequences.
//' @param cleavagePattern Regex pattern for protease cleavage sites. The default is the pattern for trypsin.
//' @return IntegerVector of missed cleavage counts.
//'
//' @examples
//' nMissed(c("PEPTIDER", "PEPKTIDER", "PEPKTIDERKK"))
//'
// [[Rcpp::export]]
Rcpp::IntegerVector nMissed(const Rcpp::CharacterVector& sequences,
                            std::string cleavagePattern = "([RK])(?=[^P])")
{
    std::regex pattern(cleavagePattern);
    size_t len = sequences.size();
    Rcpp::IntegerVector ret(len);

    for(size_t i = 0; i < len; i++) {
        std::string seq = std::string(sequences[i]);
        auto begin = std::sregex_iterator(seq.begin(), seq.end(), pattern);
        auto end = std::sregex_iterator();
        ret[i] = std::distance(begin, end);
    }

    return ret;
}

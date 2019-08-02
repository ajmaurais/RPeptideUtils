
#include <Rcpp.h>
#include <stdexcept>
#include <vector>
#include <set>
#include <map>
#include <algorithm>

#include <fastaFile.hpp>
#include <molecularFormula.hpp>
#include <peptideUtils.hpp>
#include <utils.hpp>

//!Return data files included in an R package
std::string _getPackageData(std::string filename,
														std::string packageName = "peptideUtils")
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

//' Get protein sequences for a vector of uniprot IDs
//' 
//' @title Get protein IDs in fasta file for a vector of Uniprot IDs
//' @param ids CharacterVector of uniprot IDs
//' @param fastaPath path to fasta formated file to look up protein sequences. 
//' @return CharacterVector of protein sequences in same order as ids.
//' 
//' @examples
//' #By default the fasta file included in the package containing human protein sequences is used.
//' getSquences(c("A0MZ66", "A6NMY6", "O00213", "O00213"))
//' 
//' #A fasta file can also be manually spedified.
//' fasta_path <- system.file('extdata/Human_uniprot-reviewed_20171020.fasta', package = 'peptideUtils')
//' getSquences(c("A0MZ66", "A6NMY6", "O00213", "O00213"), fasta_path)
//' 
// [[Rcpp::export]]
Rcpp::CharacterVector getSquences(const Rcpp::CharacterVector& ids, std::string fastaPath = "")
{
	std::string _fastaPath = fastaPath.empty() ? _getPackageData("extdata/Human_uniprot-reviewed_20171020.fasta") : fastaPath;
	
	Rcpp::CharacterVector ret;
	
	utils::FastaFile fasta(_fastaPath);
	if(!fasta.read()) throw std::runtime_error("Could not read fasta file!");
	
	size_t len = ids.size();
	for(size_t i = 0; i < len; i++){
		ret.push_back(fasta.getSequence(std::string(ids[i])));
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
//' @param fastaPath path to fasta formated file to look up protein sequences
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
																					const Rcpp::CharacterVector& peptideSeq,
																					std::string fastaPath = "",
																					std::string modSep = "|")
{
	std::string _fastaPath = fastaPath.empty() ? _getPackageData("extdata/Human_uniprot-reviewed_20171020.fasta") : fastaPath;
	
	size_t len = ids.size();
	if(len != peptideSeq.size())
		throw std::runtime_error("ids.size() != peptideSeq.size()");
	
	//init FastaFile
	utils::FastaFile fasta(_fastaPath);
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

//' Combined concated mods from multiple peptides into a single string.
//' 
//' @title Combined mods from multiple peptides into a single string
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
//' @param atomMasses Path to atomMasses file. If blank, the default file included in the package is used.
//' @return vector of peptide masses.
//'
//' @examples
//' calcMass(c("ACLLPETVNMEEYPYDAEY", "ALCAEFK"))
//'
// [[Rcpp::export]]
Rcpp::NumericVector calcMass(const Rcpp::StringVector& sequences,
														 bool monoMass = true,
														 std::string residueAtoms = "",
														 std::string atomMasses = "")
{
	//get data file paths
	std::string atomMassesPath = atomMasses.empty() ? _getPackageData("atomMasses.txt") : atomMasses;
	std::string residueAtomsPath = residueAtoms.empty() ? _getPackageData("defaultResidueAtoms.txt") : residueAtoms;
	char avg_mono = monoMass ? 'm' : 'a';
	
	//init residues
	utils::Residues residues(residueAtomsPath, atomMassesPath);
	if(!residues.initialize()) throw std::runtime_error("Error reading required files for calcMass!");
	
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
	if(!residues.readAtomCountTable(residueAtomsPath))
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
//' @param atomMasses Path to atomMasses file. If blank, the default file included in the package is used.
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
				  bool mz_filter = true, std::string residueAtoms = "", std::string atomMasses = "",
				  double minMz = 400, double maxMz = 1800,
				  int minCharge = 1, int maxCharge = 5,
				  size_t minLen = 6, size_t maxLen = 0)
{
	//get file paths for atom mass tables
	std::string atomMassesPath;
	std::string residueAtomsPath;
	
	if(mz_filter){
		 residueAtomsPath = residueAtoms.empty() ? _getPackageData("defaultResidueAtoms.txt") : residueAtoms;
		 atomMassesPath = atomMasses.empty() ? _getPackageData("atomMasses.txt") : atomMasses;
	}

	//check args
	if(ids.size() != sequences.size())
		throw std::runtime_error("sequences and ids must be the same length!");
	
	size_t _maxLen = maxLen == 0 ? std::string::npos : maxLen;
	Rcpp::List ret;
	
	size_t len = sequences.size();
	for(size_t i = 0; i < len; i++)
	{
		std::vector<std::string> peptides_temp;
		Rcpp::CharacterVector ret_temp;
		if(mz_filter)
		{
			//init residues
			utils::Residues residues(residueAtomsPath, atomMassesPath);
			if(!residues.initialize())
				throw std::runtime_error("Error reading required files for calcMass!");
			
			residues.digest(utils::removeWhitespace(std::string(sequences[i])), peptides_temp, nMissedCleavages, false, cleavagePattern,
							minMz, maxMz, minCharge, maxCharge);
			
			std::sort(peptides_temp.begin(), peptides_temp.end(), utils::strLenCompare());
			for(auto it = peptides_temp.begin(); it != peptides_temp.end(); ++it)
			{
				size_t len_temp = it->length();
				if(len_temp >= minLen && (_maxLen == std::string::npos ? true : len_temp <= _maxLen))
					ret_temp.push_back(it->c_str());
			}
			
		}//end if mz_filter
		else{
			utils::digest(utils::removeWhitespace(std::string(sequences[i])), peptides_temp,
						  nMissedCleavages, minLen, _maxLen, cleavagePattern);
			for(auto it = peptides_temp.begin(); it != peptides_temp.end(); ++it){
				ret_temp.push_back(it->c_str());
			}
		}
		ret.push_back(ret_temp, std::string(ids[i]));
	}
	
	return ret;
}




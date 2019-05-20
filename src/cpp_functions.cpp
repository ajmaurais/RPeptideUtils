
#include <Rcpp.h>
#include <stdexcept>
#include <vector>
#include <set>

#include <fastaFile.hpp>
#include <molecularFormula.hpp>

//Characters representing dynamic modifications
const char* MOD_CHARS = "*";

//' Get protein sequences for a vector of uniprot IDs
//' 
//' @title Get protein IDs in fasta file for a vector of Uniprot IDs
//' @param fastaPath path to fasta formated file to look up protein sequences
//' @param ids CharacterVector of uniprot IDs
//' @return CharacterVector of protein sequences in same order as ids.
//'
//' 
// [[Rcpp::export]]
Rcpp::CharacterVector getSquences(std::string fastaPath, const Rcpp::CharacterVector& ids)
{
  Rcpp::CharacterVector ret;
  
  utils::FastaFile fasta(fastaPath);
  if(!fasta.read()) throw std::runtime_error("Could not read fasta file!");
  
  size_t len = ids.size();
  for(size_t i = 0; i < len; i++){
    ret.push_back(fasta.getSequence(std::string(ids[i])));
  }
  
  return ret;
}

/**
Record indices of occurrences of MOD_CHARS in seq
\param seq peptide sequence containing modifications
\param modLocs vector populated with locations of locations

\return peptide sequence with mods removed
*/
std::string _getModLocs(std::string seq, std::vector<int>& modLocs)
{
	modLocs.clear();
	std::string ret = "";
	for(size_t i = 0; i < seq.length(); i++)
	{
		for(const char* p = MOD_CHARS; *p; p++)
		{
			if(seq[i] == *p){
				seq.erase(i, 1);
				modLocs.push_back(int(ret.length() - 1));
				break;
			}
		}
		//exit loop if at end of sequence
		if(i >= seq.length())
			break;
		
		//Check that current char is letter
		if(!isalpha(seq[i]))
			throw std::runtime_error("Invalid peptide sequence: " + seq);
		
		//add new amino acid to ret
		ret.push_back(seq[i]);
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
//' 
// [[Rcpp::export]]
Rcpp::CharacterVector getModifiedResidues(std::string fastaPath,
										  const Rcpp::CharacterVector& ids,
										  const Rcpp::CharacterVector& peptideSeq,
										  std::string modSep = "|")
{
	size_t len = ids.size();
	if(len != peptideSeq.size())
		throw std::runtime_error("ids.size() != peptideSeq.size()");
	
	//init FastaFile
	utils::FastaFile fasta(fastaPath);
  	if(!fasta.read()) throw std::runtime_error("Could not read fasta file!");

	Rcpp::CharacterVector ret(len, "");
	std::vector<int> modIndex_temp;
	for(size_t i = 0; i < len; i++)
	{
		std::string this_modLocs;
		std::string seqTemp = _getModLocs(std::string(peptideSeq[i]), modIndex_temp);
		for(auto it = modIndex_temp.begin(); it != modIndex_temp.end(); ++it){
			utils::addChar(fasta.getModifiedResidue(std::string(ids[i]), seqTemp, *it), this_modLocs, modSep);
		}
		ret[i] = this_modLocs;
	}
	return ret;
}

//' Combined concated mods from multiple peptides into a single string
//' 
//' @title Combined mods from multiple peptides into a single string
//' @param mods Modifications to combine
//' @param sep delimiter separating modifications
//' @return Modifications combined into a single string
//'
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

//' Calculate peptide monoisotopic or average masses.
//' 
//' @title Calculate peptide masses
//' @param sequences Peptide sequences
//' @param monoMass Should monoisotopic mass be calculated. If false, average mass is calculated.
//' @param residueMasses Path to residueAtoms file. If blank, the default file included in the package is used. 
//' @param atomMasses Path to atomMasses file. If blank, the default file included in the package is used.
//' @return vector of peptide masses.
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
//' @param residueMasses Path to residueAtoms file. If blank, the default file included in the package is used. 
//' @param atomMasses Path to atomMasses file. If blank, the default file included in the package is used.
//' @return vector of peptide formulas.
//'
// [[Rcpp::export]]
Rcpp::StringVector calcFormula(const Rcpp::StringVector& sequences, 
                               bool subscripts = false,
                               std::string residueAtoms = "",
                               std::string atomMasses = "")
{
  //get data file paths
  std::string atomMassesPath = atomMasses.empty() ? _getPackageData("atomMasses.txt") : atomMasses;
  std::string residueAtomsPath = residueAtoms.empty() ? _getPackageData("defaultResidueAtoms.txt") : residueAtoms;

  //init residues
  utils::Residues residues(residueAtomsPath, atomMassesPath);
  if(!residues.initialize()) throw std::runtime_error("Error reading required files for calcFormula!");
  
  size_t len = sequences.size();
  Rcpp::StringVector ret(len);
  for(size_t i = 0; i < len; i++){
    ret[i] = residues.calcFormula(std::string(sequences[i]), subscripts);
  }
  
  return ret;
}




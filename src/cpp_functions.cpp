
#include <Rcpp.h>
#include <stdexcept>
#include <vector>

#include <fastaFile.hpp>

//Characters representing dynamic modifications
const char* MOD_CHARS = "*";

//' Get protein sequences for a vector of uniprit IDs
//' 
//' @title Get protein IDs in fasta file for a vector of Uniprot IDs
//' @param fastaPath path to fasta formated file to look up protein sequences
//' @param ids CharacterVector of uniprot IDs
//' @return CharacterVector of protein sequences in same order as ids.
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

//' Get locations of modified residues in parent protein
//'
//' @title Get locations of modified residues in parent protein
//' @param fastaPath path to fasta formated file to look up protein sequences
//' @param ids CharacterVector of Uniprot IDs
//' @param peptideSeq list of peptide sequences containing modifications
//' @param modSep delimiter for multiple modifications
//' @return CharacterVector containing locations of modifications in protein sequence
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



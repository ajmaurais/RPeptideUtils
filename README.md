# peptideUtils
Wraps useful c++ functions for dealing with peptide data in R.

## Installation
`peptideUtils` depends on the [utils](https://github.com/ajmaurais/utils) c++ library which is included as a submodule in this repository. To clone this repository and the `utils` submodule, run:
```bash
git clone --recurse-submodules https://github.com/ajmaurais/peptideUtils
```
Next, build the utils library:
```bash
cd peptideUtils/utils && make
```
Finally, build and install the `peptideUtils` package.
```bash
cd ..
Rscript -e "install.packages('.', repo = NULL, type = 'source')"
```

## Features

### Peptide formula and mass calculations

```R
> calcFormula(c("ACLLPETVNMEEYPYDAEY", "ALCAEFK"), subscripts = TRUE)
[1] "C₁₀₂H₁₄₇N₂₁O₃₆S₂" "C₃₇H₅₉N₉O₁₁S"

> calcMass(c("ACLLPETVNMEEYPYDAEY", "ALCAEFK"))
[1] 2305.9764  837.4057
```

### Work with peptide modifications

```R
# Look up locations of modified residues in parent protein
> getModifiedResidues(c("Q00839", "Q9HCS7", "Q7L014"), c("APQC*LGK", "FADMEC*K", "GAEIIVC*TPGR"))
[1] "C562" "C676" "C501"

# Concatenate a vector of modified residues into a single string
> combineMods(c('C157', 'C157|C125', 'C50', 'C125'))
[1] "C125|C157|C50"
```

### Look up protein sequences by Uniprot ID
```R
> getSquences(c("A0MZ66", "A6NMY6"))
[1] "MNSSDEEKQLQLITSLK..."
[2] "MSTVHEILCKLSLEGDH..."  
```

### Convert between 1 and 3 letter amino acid codes

```R
> oneLetterToThree(c("AC*LLPETVNMEEYPYDAEY", "ALCAEFK", "AQUPIVER"))
[1] "AlaCys*LeuLeuProGluThrValAsnMetGluGluTyrProTyrAspAlaGluTyr"
[2] "AlaLeuCysAlaGluPheLys"                                  
[3] "AlaGlnSecProIleValGluArg"

> threeLetterToOne(c("Ala-Cys*-Leu-Leu-Pro", "Ala-Leu-Cys-Ala", "Ala-Gln-Sec-Ile"), sep_in="-")
[1] "AC*LLP" "ALCA"   "AQUI"
```

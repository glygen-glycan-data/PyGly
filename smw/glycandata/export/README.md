# Documentation for GlyGen Glycan Data export files

## General

Header field GlyTouCanAccession consistently refers to GlyTouCan accession 	

### File-extensions

- file.tsv  Tab-separated values, with headers/field names in the first line
- file.txt  Other text-based format
- file.zip  Zip-file of gtcacc.txt for each GlyTouCan accession (typically sequence formats). 
  
## monocounts.tsv

Triples representing monosaccharide counts, computed by the Glycan.iupac_composition method in PyGly. 

Fields
- accession	- GlyTouCan accession
- monosaccharide - IUPAC symbol for monosaccharide (also Xxx, X, \*)
- count - integer count for monosaccharide symbol. Counts for monosaccharides contained in repeated units are marked with a "+"

Notes
 - IUPAC symbols with sterochemistry (e.g. Glc, Gal, Man) have their counts aggregated to their subsuming IUPAC symbol (e.g. Hex). 
- Xxx represents any monosaccharide not otherwise captured by an IUPAC symbol. 
- X represents any (floating) substituent not otherwise captured by an IUPAC symbol. 
- \* represents the total number of monosaccharides
- \+aldi monosaccharides (e.g. Glc+aldi) are reported, but also as their components.
- Substituents S, P, Me are split from their residues and counted separately.
- Other floating substituents are counted as X
- Non-floating substituents are left attached to their residue, resulting in Xxx (except for GlcNAc, etc.)


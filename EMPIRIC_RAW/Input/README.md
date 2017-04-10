
##### Oligo sequences for creating EMPIRIC libraries: Primers and cassettes.xlsx

##### Input files for filtering deep sequencing data are grouped by ortholog.

##### Summary of fitness values used for analysis: db-fitness.csv
Fitness values are normalized to the average stop of each library. Raw sequencing files (fastq) available by email. 


#### ARCHIVE:  IGPS-paper-rawdata.tar and IGPS-paper-replicate.tar

IGPS-paper-rawdata.tar: Compressed file with raw data files (counts AA/time and codon/time) including all accessory files

IGPS-paper-replicate.tar: Compressed file with raw data files (counts AA/time and codon/time) for biological replicates of SsIGPS B3 and B4 libraries
  -  Used to demostrate reproducibility of fitness results.  

*************************************************************** 
RAW MUTANT COUNTS FOR IGPS PROTEINS USED IN PUBLICATION: 

Correlation of fitness landscapes from three orthologous TIM barrels originates from sequence and structure constraints

Y.H. Chan, S.V. Venev, K.B. Zeldovich, C.R. Matthews

Nature Communications 8:14614 (2017) DOI: 10.1038/ncomms14614

*************************************************************** 

For each protein, Ss, Tm, Tt there is a folder with the following files:

PDB file with the structures
Some residues have been changed compared to the original PDB from Protein databank, to reflect the actual plasmid sequence
Ss: R43S
Tm: deleted residues 1-31, removed chain B, mutated C102S
Tt: deleted residues 1-34, removed chain B, no mutations.

Fasta file with the nucleotide sequence of the gene as on the plasmid

**Mapping tables for libraries 1 through 8. 
Fields: library name, library position (0-9), WT amino acid in the gene, WT codon in the gene, PDB residue number, PDB amino acid
The PDB amino acid (from the PDB provided) matches the WT amino acid in the plasmid.
Library position numbering goes 0-9 within each library. PDB residue number is as in PDB file.

**Table with read counts for each mutant: outtable-aa-bytime.txt
Fields: library name, mutant, timepoint (hours), count
Library names are b1 through b8 for mutants at different timepoints, and b1-wt through b8-wt for sequencing of the WT plasmid to assess sequencing errors. For b1-wt through b8-wt, there is only one timepoint, zero.

** Lists of MmeI sites for each library. These mutations have at least one mutant codon creating a nuclease cut site, and so the mutant counts are either depleted or less reliable.

** Blacklist –  for Tm and Tt proteins, list of libraries/timepoint combinations with consistently have very low sequencing quality and/or depth, which we removed from analysis




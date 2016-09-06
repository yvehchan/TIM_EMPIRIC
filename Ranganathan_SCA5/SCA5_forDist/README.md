- contents of this folder must be placed inside the folder with the SCA toolbox installation (downloaded from http://systems.swmed.edu/rr_lab/sca.html), such that ./sca5 is on the same level as IGPS.m, and the IGPS-input files.

- IGPS-input files include, 1IGS.pdb (Ss IGPS protein structure), and ncbi-nr-537-alignmen_cleaned.free (special format MSA file)

- IGPS.m is the modified script with the IGPS-SCA analysis. It is based on the Tutorial_pdz_v5.m script that is distributed with the SCA v5.0

- preprocess_fasta_MSA.py script does some simple MSA preprocessing: given an alignment file in fasta format, the script will eliminate columns with >80% of gaps, saving _cleaned.fasta MSA variant, and then it’ll also save such “cleaned“ alignment in a free format used by IGPS.m script.

- generated pml file can be viewed in PyMol.

# TIM_EMPIRIC

Assorted constellation of scripts and procedure for processing and analysis of the EMPIRIC experiment performed on three IGPS orthologs. 


The following is a list of the folders' descriptions:

- ./Inputs - this folder contains preprocessed EMPIRIC fitness readouts and structural/functional description of the mutated positions.

- ./numbers_for_paper - this folder contains some of the processing scripts to store the results in.

- ./PFam contains everything related IGPS family alignment, obtained from PFam

- ./Ranganathan_SCA5 is a folder created by unpacking an archive of SCA v5.0 toolbox for MATLAB, as downloaded from http://systems.swmed.edu/rr_lab/sca.html . Modified scripts and input files used in our analysis are provided, while the original contents of SCA-distribution are removed to comply with "UT SOUTHWESTERN Internal Use License Agreement". Simply aquire a copy of SCA toolbox, place modified scripts/inputs into internal folder and run the analysis.

- ./SPalignResults contains everything related structural alignment of ~70 representative TIM-barrel folds. Put SPalignNS executable with the required pdb files inside this folder and run provided scripts to obtain the Multiple Structural Sequence Alignment (MSSA).

- Analysis_Full.py is a master script for performing PCA on the EMPIRIC data, IGPS-MSA and TIM-MSA. It relies on 2 additional modules: EMPIRIC_analysis_EXP.py for EMPIRIC PCA, and MSA_analysis_EXP.py for MSAs PCA. Master script does generate some output files with principal components and such. 

- EMPIRIC_analysis_EXP.py is a Python module needed for Analysis_Full.py

- MSA_analysis_EXP.py is a Python module needed for Analysis_Full.py

- FINAL_ANALYSIS_SCRIPT_includezero_groupedAA.py is a script that generates distributions of Pearson correlation coefficient between amino acids substitution responses for pairs of residues grouped according to different criteria. Results are stored in ./numbers_for_paper .

- FINAL_ANALYSIS_SCRIPT_includezero_similar_res.py is a script that generates distributions of Pearson correlation coefficient between amino acids substitution responses for pairs of residues grouped according to different criteria. (IT IS PROBABLY SIMPLY AND OLDER VERSION OF _groupedAA SCRIPT, SO CHECK IT ...) Results are stored in ./numbers_for_paper .

- PVALUE_ratio_ANALYSIS_SCRIPT_yc.py is a script that does the permutation tests for comparison of fitness distributions in different groups of residues.  Results are stored in ./
















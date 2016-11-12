# TIM_EMPIRIC postprocessing

Assorted constellation of scripts and procedures for processing and analysis of the EMPIRIC experiment performed on three IGPS orthologs. 


The following is a list of the folders' contents:

- ./Inputs 
	* This folder contains processed EMPIRIC fitness readouts and structural/functional description of the mutated positions.

- ./MSA_inputs_all
	* This folder contains all the multiple sequence alignments (MSA) used for analysis. 

- ./numbers_for_paper 
	* This folder contains the output for some of the scripts.

- ./PFam 
	* This folder contains files obtained from PFam related to the IGPS family alignment

- ./SPalignResults
    * This folder contains files related to the structural alignment of ~70 representative TIM-barrel folds. 
	* To obtain the Multiple Structural Sequence Alignment (MSSA), put SPalignNS executable with the required pdb files inside this folder and run provided scripts.
	
- ./Ranganathan_SCA5 
	* This folder is an archive of the SCA v5.0 toolbox for MATLAB, as downloaded from http://systems.swmed.edu/rr_lab/sca.html 
	* Modified scripts and input files used in our analysis are provided, while the original contents of SCA-distribution are removed to comply with "UT SOUTHWESTERN Internal Use License Agreement". 
	* To run our analysis, acquire a copy of SCA toolbox and place modified scripts/inputs into the resulting folder.

- ./Plots 
	* DFE and Correlation_analyses
		* R scripts to generate representative graphs of dataset 
	* Structural analyses
		* Pymol commands to generate structural plots

Scripts: 
- Analysis_Full.py 
	* Master script for performing PCA on the EMPIRIC data, IGPS-MSA and TIM-MSA. 
	* Two dependent modules are required to run this script: EMPIRIC_analysis_EXP.py for EMPIRIC PCA, and MSA_analysis_EXP.py for MSAs PCA. 
	* Multiple outputs include principal components and biplots to relate features to principal components. 

- EMPIRIC_analysis_EXP.py 
	* Python module needed for Analysis_Full.py

- MSA_analysis_EXP.py 
	* Python module needed for Analysis_Full.py

- FINAL_ANALYSIS_SCRIPT_includezero_groupedAA.py 
	* Python script that generates distributions of Pearson correlation coefficient between amino acids substitution responses for pairs of residues grouped according to different criteria. 
	* Results are stored in ./numbers_for_paper 

- PVALUE_ratio_ANALYSIS_SCRIPT.py 
	* Python script that does the permutation tests for comparison of fitness distributions in different groups of residues.  Results are stored in ./




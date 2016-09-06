# EMPIRIC data processing

This folder contains a pipeline for processing fastq files. Short sequence reads are filtered for high quality and tabulated. Fitness values are calculated for different nucleotides and amino acid substitutions based on the sequence counts.


- analyze-empiric.pl is a Perl wrapper-script with multiple dependecies from ./src folder, including compiled "processfasta" C-program.

- input files to included later (FASTQs, PDBs, primer and barcode sequences, etc.)

- Original folder organization is kept to minimize adoption and reading issues.






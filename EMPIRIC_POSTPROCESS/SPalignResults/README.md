- functioning executable of the SPalignNS must be inside of this directory. SPalignNS vs SPalign does not seem to make much of a difference, because “-NS” flag determines if NS-algorithm is turned on or not (check that more carefully).

- WtIGPS-pos provides a list of positions analyzed in the libraries. These positions are extracted from pairwise structural alignments.

- run_pdb launches SPalign ~70 times for each of the templates Ss,Tm,Tt.

- parse_spalign_aln parses pairwise alignments and extracts library positions. A combined alignment is produced at the end.

- Xx_tim_aln.fasta are the resulting combined alignments.

- corresponding files must be stored in a ./pdbA directory.

- templates are:
Tm: pdb1i4n_A.ent
Tt: pdb1vc4_A.ent
Ss: pdb2c3z_A.ent

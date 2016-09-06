- functioning executable of the SPalignNS must be inside of this directory. SPalignNS vs SPalign does not seem to make much of a difference, because “-NS” flag determines if NS-algorithm is turned on or not (check that more carefully).

- WtIGPS-pos, provides a list of positions analyzed in the library, they are to be extracted from pairwise structural alignments.

- run_pdb launches SPalign ~70 times for each of the templates Ss,Tm,Tt.

- parse_spalign_aln parses pairwise alignments and extracts library positions, combined alignment is produced at the end.

- Xx_tim_aln.fasta are the resulting combined alignments.


- list of 74 representative TIM-barrel structure is provided:
pdb1a0c_A.ent  pdb1aq0_A.ent  pdb1bd0_A.ent  pdb1byb_A.ent  pdb1dl3_A.ent  pdb1eex_A.ent  pdb1gow_A.ent  pdb1huv_A.ent  pdb1jcl_A.ent  pdb1km0_A.ent  pdb1onr_A.ent  pdb1qnr_A.ent  pdb1tml_A.ent  pdb2c3z_A.ent  pdb4ubp_A.ent
pdb1a3x_A.ent  pdb1b3o_A.ent  pdb1bf6_A.ent  pdb1ceo_A.ent  pdb1dos_A.ent  pdb1eom_A.ent  pdb1gvf_A.ent  pdb1i1w_A.ent  pdb1juk_A.ent  pdb1l2q_A.ent  pdb1pdy_A.ent  pdb1qo2_A.ent  pdb1ttp_A.ent  pdb2dor_A.ent  pdb5ptd_A.ent
pdb1ad1_A.ent  pdb1b4k_A.ent  pdb1bpl_A.ent  pdb1cnv_A.ent  pdb1e70_M.ent  pdb1f74_A.ent  pdb1gw1_A.ent  pdb1i4n_A.ent  pdb1k4g_A.ent  pdb1luc_A.ent  pdb1pym_A.ent  pdb1qtw_A.ent  pdb1uro_A.ent  pdb2ebn_A.ent  pdb7a3h_A.ent
pdb1ado_A.ent  pdb1b54_A.ent  pdb1bqc_A.ent  pdb1d3g_A.ent  pdb1edg_A.ent  pdb1f8m_A.ent  pdb1h7w_A.ent  pdb1itq_A.ent  pdb1kbl_A.ent  pdb1lwh_A.ent  pdb1qap_A.ent  pdb1req_A.ent  pdb1vc4_A.ent  pdb2tmd_A.ent  pdb8ruc_A.ent
pdb1ads_A.ent  pdb1b5t_A.ent  pdb1btm_A.ent  pdb1d8c_A.ent  pdb1edq_A.ent  pdb1gg0_A.ent  pdb1hg3_A.ent  pdb1j5s_A.ent  pdb1kd0_A.ent  pdb1nar_A.ent  pdb1qat_A.ent  pdb1smd_A.ent  pdb1vrx_A.ent  pdb2tps_A.ent

- corresponding files must be stored in a ./pdbA directory.

- templates are:
Tm: pdb1i4n_A.ent
Tt: pdb1vc4_A.ent
Ss: pdb2c3z_A.ent

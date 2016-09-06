import re
import sys
import pandas as pd
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
import subprocess as sub


# tm_pdb = ('Tm',"pdb1i4n_A.ent")
# tt_pdb = ('Tt',"pdb1vc4_A.ent")
# ss_pdb = ('Ss',"pdb2c3z_A.ent")

if   sys.argv[1] == 'Tm':
    template_name, current_template, template_offset = 'Tm',"pdb1i4n_A.ent", -1
elif sys.argv[1] == 'Tt':
    template_name, current_template, template_offset = 'Tt',"pdb1vc4_A.ent", 0
elif sys.argv[1] == 'Ss':
    template_name, current_template, template_offset = 'Ss',"pdb2c3z_A.ent", -28
else:
    print "blah! wrong input: Tm Tt or Ss are acceptable only!"
    print "nohup python run_pdb.py Tm(Tt,Ss) &"
    sys.exit(1)

fnames = [s.strip() for s in sub.check_output('ls ./ | grep "_ent_%s.aln"'%template_name, shell=True).strip().split('\n')]
# # # # read all the fnames ...
# # fnames = [f.strip() for f in sys.stdin.readlines()]
# template_name = 'Tm'
# template_offset = -1


def parse_sp_aln(sp_aln_fname):
    fp = open(sp_aln_fname,'r')
    content = fp.readlines()
    fp.close()
    #
    # figure out match offset ...
    idx = 0
    for char in content[12]:
        if char in ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']:
            break
        idx += 1
    match_offset = idx
    #
    tmp_start,tmp_seq,tmp_stop = re.split( '\s+', content[12].strip() )
    aln_match = content[13][match_offset:]
    aln_start,aln_seq,aln_stop = re.split( '\s+', content[14].strip() )
    return aln_match,(tmp_start,tmp_seq,tmp_stop),(aln_start,aln_seq,aln_stop)


# make a gapped--ungapped map ...
def get_aln_tmp_map(full_aln_info,offset,lib_pos):
    # unpack stuff ...
    match,tmp_info,aln_info = full_aln_info
    start,tmp_seq,stop = tmp_info
    _,aln_seq,_ = aln_info
    # correct real numbering to PDB-trimmed numbering ...
    lib_pos_pdb = [pos+offset for pos in lib_pos]
    # map positions in PDB to gapped alignment ...
    ungapped_idx = int(start) # fix this ...
    tmp_aln_map = []
    for idx,aa in enumerate(tmp_seq):
        if (aa != '-'):
            ungapped_idx += 1
            tmp_aln_map.append( (ungapped_idx,idx) )
    # turn map into dictionary ...
    the_map = dict(tmp_aln_map)
    #
    seq_to_comape = []
    matches = []
    s2 = []
    for pos in lib_pos_pdb:
        if pos in the_map.keys():
            pos_in_aln = the_map[pos]
            seq_to_comape.append(tmp_seq[pos_in_aln])
            matches.append(match[pos_in_aln])
            s2.append(aln_seq[pos_in_aln])
        else:
            seq_to_comape.append('-')
            matches.append(' ')
            s2.append('-')

    ######
    # print
    print ''.join(seq_to_comape)
    # print ''.join(matches)
    print ''.join(s2)
    print ''.join(matches)
    print
    return ''.join(s2)


lib_dat = pd.read_csv("WtIGPS-pos.csv")
tmp_dat = lib_dat.groupby('organism').get_group(template_name)
lib_pos = list(tmp_dat['pos'])
lib_aas = list(tmp_dat['wtaa'])

# print ''.join(lib_aas)

final_lib_aln = []
for fname in fnames:
    print fname
    print ''.join(lib_aas)
    sss = get_aln_tmp_map( parse_sp_aln(fname), template_offset, lib_pos=lib_pos)
    final_lib_aln.append(sss)


def fname_to_id(fname):
    fname_major = fname.split('_')[0]
    if   fname_major == 'pdb1i4n':
        return 'Tm'
    elif fname_major == 'pdb1vc4':
        return 'Tt'
    elif fname_major == 'pdb2c3z':
        return 'Ss'
    else:
        return fname.strip('.aln')    

# wrap generated alignments into list of SeqRec...
seq_recs_to_output = [ SeqRecord.SeqRecord( Seq.Seq(seq), id=fname_to_id(fname), description='') for seq,fname in zip(final_lib_aln,fnames) ]
# output ...
SeqIO.write( seq_recs_to_output, "%s_tim_aln.fasta"%template_name, 'fasta' )


# # ################## Alignment Report ##################
# # pdb1i4n_A.ent pdb1a0c_A.ent: Length= 251 437
# # Pfold= 54.5 %; SPe/SPa/SPb= 0.531 0.497 0.620 ;Effective_Length: 313
# # RMSD/Nali= 3.40 / 192 ;GDT= 0.406 ;TMscore(a,b,c)= 0.469 0.612 0.381 SEQID= 8.0%

# # Rotation Matrix:
# #   55.66567 -0.09830  0.37672 -0.92110
# #   84.25718 -0.73650 -0.65001 -0.18724
# #   45.87282 -0.66926  0.65998  0.34135

# # Alignment: pdb1i4n_A.ent pdb1a0c_A.ent        Structure Overlap: 54.18 %
# # (':' denotes the residue pairs of distance <= 4A, and '.' denotes <=8A)
# # RRLWEIVEAKKKDILEIDGENLIVQRRNHRFLEVLSGKE----RVKIIAEFK---KASPSA---GDINAD-A-----------SL-EDFIRMYDE--LADAISILTEKHYFKG-D------------P-AFVR-AARNLTC----RPILAK-DF---------YI----------DTVQVKLA-SSV-GA-D-AILIIAR-I-------------L--TAEQIKEIYEAAEEL-GMD-SLVEVH-----------S-REDLEKVFS-VI-R-PKIIGIN-TRDLDTFEIKK--NVLWELLPLV--PDDTVVVAE-SG----------IKDPRE-LKD-LRGKV-----N---AVLV-GTS-IMK---AENPRR-FLEE-MRA--W--SE----------------------------------------------------------------------
# #                                       .    :::::::::   ::. .:   .....: .           .: ::::::: :  ::::::::::     : :            : .. : ::::::.    ::.::: ::         ..          :::::::: :.. :: : .::::.. :             .  ::::::::::::::: :.: ::::::           : ::::::::: :: : : ::::: ::      .:.  :.::::::::  :: :.:::: .:          .::..: .:: :::::     :   :::: :::   .   .:.: . :. . : .  :  ..
# # ----NKYFENVSKIKYEGPKSNNPYSFKFYNPEEVIDGKTMEEHLRFSIAYWHTFTADGTDQFGKATMQRPWNHYTDPMDIAKARVEAAFEFF-DKINAPYFCFHDR-----DIAPEGDTLRETNKNLDTI-VAMIKDYLKTSKTKVLWGTANLFSNPRFVHGASTSCNADVFAYSAAQVKKALEITKELGGENYVFWGGREGYETLLNTDMEFELDNFARFLHMAVDYAKEIGFEGQFLIEPKPKEPTKHQYDFDVANVLAFLRKYDLDKYF-KVNIEANH------ATLAFHDFQHELRYARING-VLGSIDANTGDMLLGWDTDQFPTDIRMTTLAMYEVIKMGGFDKGGLNFDAKVRRASFEPEDLF-LGHI-AGM-DAFAKGFKVAYKLVKDRVFDKFIEERYASYKDGIGADIVSGKADFRSLEKYALERSQIVNKSGRQELLESILNQYLFA



# 4    EIVEAKKKDILEIDGENLIVQRRNHRFLEVLSGKE----RVKIIAEFK---KASPSA---GDINAD-A-----------SL-EDFIRMYDE--LADAISILTEKHYFKG-D------------P-AFVR-AARNLTC----RPILAK-DF---------YI----------DTVQVKLA-SSV-GA-D-AILIIAR-I-------------L--TAEQIKEIYEAAEEL-GMD-SLVEVH-----------S-REDLEKVFS-VI-R-PKIIGIN-TRDLDTFEIKK--NVLWELLPLV--PDDTVVVAE-SG----------IKDPRE-LKD-LRGKV-----N---AVLV-GTS-IMK---AENPRR-FLEE-MRA--W--SE 250
#                                        .    :::::::::   ::. .:   .....: .           .: ::::::: :  ::::::::::     : :            : .. : ::::::.    ::.::: ::         ..          :::::::: :.. :: : .::::.. :             .  ::::::::::::::: :.: ::::::           : ::::::::: :: : : ::::: ::      .:.  :.::::::::  :: :.:::: .:          .::..: .:: :::::     :   :::: :::   .   .:.: . :. . : .  :  ..
# 0    NKYFENVSKIKYEGPKSNNPYSFKFYNPEEVIDGKTMEEHLRFSIAYWHTFTADGTDQFGKATMQRPWNHYTDPMDIAKARVEAAFEFF-DKINAPYFCFHDR-----DIAPEGDTLRETNKNLDTI-VAMIKDYLKTSKTKVLWGTANLFSNPRFVHGASTSCNADVFAYSAAQVKKALEITKELGGENYVFWGGREGYETLLNTDMEFELDNFARFLHMAVDYAKEIGFEGQFLIEPKPKEPTKHQYDFDVANVLAFLRKYDLDKYF-KVNIEANH------ATLAFHDFQHELRYARING-VLGSIDANTGDMLLGWDTDQFPTDIRMTTLAMYEVIKMGGFDKGGLNFDAKVRRASFEPEDLF-LGHI-AGM-DAFAKGFK 366

# fnames = ['pdb1a0c_A.ent.aln', 'pdb1a3x_A.ent.aln', 'pdb1ad1_A.ent.aln', 'pdb1ado_A.ent.aln', 'pdb1ads_A.ent.aln',
# 'pdb1aq0_A.ent.aln', 'pdb1b3o_A.ent.aln', 'pdb1b4k_A.ent.aln', 'pdb1b54_A.ent.aln', 'pdb1b5t_A.ent.aln', 'pdb1bd0_A.ent.aln',
# 'pdb1bf6_A.ent.aln', 'pdb1bpl_A.ent.aln', 'pdb1bqc_A.ent.aln', 'pdb1btm_A.ent.aln', 'pdb1byb_A.ent.aln', 'pdb1ceo_A.ent.aln', 'pdb1cnv_A.ent.aln',
# 'pdb1d3g_A.ent.aln', 'pdb1d8c_A.ent.aln', 'pdb1dl3_A.ent.aln', 'pdb1dos_A.ent.aln', 'pdb1e70_M.ent.aln', 'pdb1edg_A.ent.aln',
# 'pdb1edq_A.ent.aln', 'pdb1eex_A.ent.aln', 'pdb1eom_A.ent.aln', 'pdb1f74_A.ent.aln', 'pdb1f8m_A.ent.aln', 'pdb1gg0_A.ent.aln',
# 'pdb1gow_A.ent.aln', 'pdb1gvf_A.ent.aln', 'pdb1gw1_A.ent.aln', 'pdb1h7w_A.ent.aln', 'pdb1hg3_A.ent.aln', 'pdb1huv_A.ent.aln',
# 'pdb1i1w_A.ent.aln', 'pdb1i4n_A.ent.aln', 'pdb1itq_A.ent.aln', 'pdb1j5s_A.ent.aln', 'pdb1jcl_A.ent.aln', 'pdb1juk_A.ent.aln',
# 'pdb1k4g_A.ent.aln', 'pdb1kbl_A.ent.aln', 'pdb1kd0_A.ent.aln', 'pdb1km0_A.ent.aln', 'pdb1l2q_A.ent.aln', 'pdb1luc_A.ent.aln',
# 'pdb1lwh_A.ent.aln', 'pdb1nar_A.ent.aln', 'pdb1onr_A.ent.aln', 'pdb1pdy_A.ent.aln', 'pdb1pym_A.ent.aln', 'pdb1qap_A.ent.aln',
# 'pdb1qat_A.ent.aln', 'pdb1qnr_A.ent.aln', 'pdb1qo2_A.ent.aln', 'pdb1qtw_A.ent.aln', 'pdb1req_A.ent.aln', 'pdb1smd_A.ent.aln',
# 'pdb1tml_A.ent.aln', 'pdb1ttp_A.ent.aln', 'pdb1uro_A.ent.aln', 'pdb1vc4_A.ent.aln', 'pdb1vrx_A.ent.aln', 'pdb2c3z_A.ent.aln',
# 'pdb2dor_A.ent.aln', 'pdb2ebn_A.ent.aln', 'pdb2tmd_A.ent.aln', 'pdb2tps_A.ent.aln', 'pdb4ubp_A.ent.aln', 'pdb5ptd_A.ent.aln',
# 'pdb7a3h_A.ent.aln', 'pdb8ruc_A.ent.aln']






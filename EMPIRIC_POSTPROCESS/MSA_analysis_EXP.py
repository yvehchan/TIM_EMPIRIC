
from Bio import AlignIO
import pandas as pd
import copy
import matplotlib.pyplot as plt
from Bio import Data
# aacids = list(Data.IUPACData.protein_letters)
from Bio import SeqUtils
from Bio.SeqUtils import ProtParamData
KD = ProtParamData.kd
#
import numpy as np
from scipy import stats as st
from Bio.Alphabet import generic_protein
# from Bio.Align import AlignInfo
from inspect import isfunction
#
aacids = sorted(list(SeqUtils.IUPAC.protein.letters))


# # EMPIRIC_raw_data_fname = "db-fitness.csv"
# # most important reference file with all the libraries information ...
# EMPIRIC_features_fname = "features-original.csv"
# # 20 amino acids ...aacids = list(SeqUtils.IUPAC.protein.letters)
# aln_fname = "alignment432_clust.fasta_ali"


def readAlnIdentifyTemplates(fname,template_ids=['SsIGPS','TmIGPS','TtIGPS']):
    """function read the alignment from file and finds template sequences in it.
    Returns alignment itself and dict with template indexes as values."""
    aln = AlignIO.read(fname,'fasta',alphabet=generic_protein)
    # check if templates are duplicated in the aln ...
    #
    #
    tmp_aln = {}
    for tmp in template_ids:
        # we are asserting that there is just one instance of each tmp in aln
        tmp_aln[tmp], = [idx for idx,seqrec in enumerate(aln) if seqrec.id==tmp]
        print "%s is located at following position in alignment: "%tmp, tmp_aln[tmp]
    #
    return aln,tmp_aln




def getAlnForTemplateLib(aln,tmpid,tmplib_pos,save_aln=False,sub_aln_fname=''):
    """Function takes global aln, organism name (template),
    and positions in library corresponding to the template."""
    ##################################
    # make a gapped--ungapped map ...
    ##################################
    def get_gap_ungap_map(gapped_seq):
        map_list = []
        ungapped_counter_one = 1
        for gapped_counter_zero,aa in enumerate(gapped_seq):
            # count non-gap AAcids ...
            if (aa != '-'):
                map_list.append((ungapped_counter_one,gapped_counter_zero))
                ungapped_counter_one += 1
        # turn map into dictionary ...
        the_map = dict(map_list)
        return the_map
    ################################
    # get sub-alignment from aln, given a set of positions of interest ...
    def get_subalignment(positions,aln):
        positions = list(positions)
        pos1 = positions[0]
        sub_aln = aln[:,pos1:pos1+1]
        for pos in positions[1:]:
            sub_aln += aln[:,pos:pos+1]
        return sub_aln
    ################################
    # create a map for THE template 
    tmp_map = get_gap_ungap_map(aln[tmpid])
    # get the positions in the alignment that correspond to the library positions ... 
    lib_pos_in_aln = [tmp_map[pos] for pos in tmplib_pos]
    # extract tmp aa just to check correctness ...
    lib_aas_in_aln = [aln[tmpid].seq[pos] for pos in lib_pos_in_aln]    
    # extract alignment columns for positions 'lib_pos_in_aln' ...
    aln_for_lib = get_subalignment(lib_pos_in_aln,aln)
    # save alignment if you wish ...
    if save_aln:
        AlignIO.write([aln_for_lib,],sub_aln_fname,'fasta')
    #################################
    return lib_pos_in_aln, lib_aas_in_aln, aln_for_lib



#################################################################################################
# make a function to get aln summary ...
#################################################################################################
def getSummary(aln):
    #
    def aa_consensus(aln_col):
        # neglect gaps altogether ...
        aln_col_filtered = aln_col.replace('-','')
        # turn into pandas Series to use its power ...
        aa_counts = pd.Series( list(aln_col_filtered) ).value_counts()
        # consensus stuff and IC ...
        aa_freq = aa_counts/float(aa_counts.sum())
        perc_cons = aa_freq.max()*100.0
        IC = sum( freq*np.log(freq) for freq in aa_freq )
        # return AA with the highest counts ...
        return aa_counts.idxmax(), perc_cons, IC
    # get alignment length, presumably ~80
    aln_length = aln.get_alignment_length()
    # gaps profiles - how many gaps are in a certain position ...
    gaps_profile = [ aln[:,pos].count('-') for pos in range(aln_length) ]
    # get consensus stuff ...
    cons_info = [ aa_consensus(aln[:,pos]) for pos in range(aln_length) ]
    # extract consensus infos ...
    consensus = [ aa for aa,perc,ic in cons_info ]
    info_cont = [ ic for aa,perc,ic in cons_info ]
    perc_cons = [ perc for aa,perc,ic in cons_info ]
    #
    return {'gaps_profile': gaps_profile,
            'consensus':    consensus,
            'ic':           info_cont,
            'perc_cons':    perc_cons}





#################################################################################################
# make a function to get AA usage matrix from Aln ...
#################################################################################################
def alnToAAmat(aln, method='counts',pos_labels=None,log_transform=True):
    # function to get aa counts from aln position ...
    def get_aa_counts(aln_col):
        # # neglect gaps altogether ...
        # aln_col_filtered = aln_col.replace('-','')
        # # # turn into pandas Series to use its power ...
        # # # value counts takes into account things like '-' and 'X', which is unwanted for our purposes ...
        # # aa_counts = pd.Series( list(aln_col_filtered) ).value_counts().sort_index()
        aa_counts = pd.Series( (str(aln_col).count(aa) for aa in aacids), index=aacids )
        return aa_counts
    # function to get aa fractions from aln position ...
    def get_aa_fractions(aln_col):
        # # neglect gaps altogether ...
        # aln_col_filtered = aln_col.replace('-','')
        # # turn into pandas Series to use its power ...
        # # # # value counts takes into account things like '-' and 'X', which is unwanted for our purposes ...
        aa_counts = pd.Series( (str(aln_col).count(aa) for aa in aacids), index=aacids )
        aa_freqs = aa_counts/aa_counts.sum()
        return aa_freqs
    #############################
    # get alignment length, presumably ~80
    aln_length = aln.get_alignment_length()
    # go through all aln positions (columns) and calculate aa counts ...
    if method=='counts':
        dat_mat = [ get_aa_counts(aln[:,pos]) for pos in range(aln_length) ]
    elif method=='fractions':
        dat_mat = [ get_aa_fractions(aln[:,pos]) for pos in range(aln_length) ]
    # the matrix with labels or without ...
    matrix = pd.DataFrame(dat_mat) if pos_labels is None else pd.DataFrame(dat_mat,index=pos_labels)
    matrix = matrix.fillna(0)
    print matrix
    # do log trans. if specified (only for counts though) ...
    if (method=='counts')and(log_transform):
        matrix = np.log(matrix+1)
    # returning ...
    return matrix



############################################################################################################################
# #add features to table
# def pivotAddFeatures(dat, dat_features):
#     """Takes raw EMPIRIC data (unrolled and filtered), and features of positions in library:
#     Returns a merged table, where each library position is characterized by features and 20-EMPIRIC fitness readings.
#     BEWARE:
#     Unique labels are 'organism-pos' and their order as in features table, which is cruicial for compatibility."""
#     #pivot EMPIRIC dataset to have mutaa(fitness readings) as columns and wt pos as rows
#     table = dat.pivot(index='organism-pos',columns='mutaa',values='fitness')
#     # add wtaa residues to table based on org-pos
#     columns_to_merge = dat[['organism-pos','wtaa']].drop_duplicates()
#     merged_table = table.merge(columns_to_merge, how='inner',left_index=True,right_on='organism-pos').reset_index(drop=True)
#     # add other features of interest to the table
#     # MUST merge the table to features to keep org-pos order as in EMPIRIC_features_fname ...
#     merged_table = dat_features.merge(merged_table,how='inner',on='organism-pos').reset_index(drop=True)
#     #add kd of wildtype
#     KD = ProtParamData.kd
#     merged_table['wtKD'] = merged_table['wtaa'].map(KD)
#     # we can add more features later on here ...
#     # ...
#     # ...
#     #print merged_table
#     return merged_table



#main function for running PCA, calls on subfunctions
def runPCA(dat,matrix_w=None):
    """ run PCA, notes:
    matrix has to be pandas df, and can contain NAs, they are ommited during covariation calculation
    and NAs are filled with 0.0 during the projection onto eigenvectors. """
    #######################################################
    def calculateVarianceExplainedSort(eig_vals, eig_vecs):
        """function zips eig val and vecs and sort them simultaneously,
        calculates fractions of explained variance and returns fracs and eigvecs
        sorted by eigvalues in descending order."""
        #####################################
        #  check if there are any negative eigvals, there should not be any presumably ...
        if (eig_vals < 0.0).any():
            print "BEWARE: There are some negative eigvals in the PCA analysis!"
            print "script proceeds, but that's something to check!"
        ###############################################
        #sort from largest to smallest eigenvalues (both eigvals and eigvectors)
        sorted_eigs = sorted(zip(np.abs(eig_vals),np.transpose(eig_vecs)),reverse=True)
        # extracted sorted vectors ...
        sorted_eigvecs = [eigvec for eigval,eigvec in sorted_eigs]
        # calculate var fractions (ordered the same way as vectors: descending)
        total_eigval = sum([eigval for eigval,eigvec in sorted_eigs])
        fracs_var_explained = [eigval/total_eigval*100.0 for eigval,eigvec in sorted_eigs]
        return (fracs_var_explained, sorted_eigvecs) 
    ################################ 
    def getProjectionMatrix(sorted_eigvecs):
        """takes eigen vectors in a sorted order and stacks them to create a projection matrix W"""
        matrix_w = np.vstack(sorted_eigvecs).transpose()
        return matrix_w
    ############################################
    def getDotProduct(origina_matrix, matrix_w):
        Y = np.asarray(origina_matrix.fillna(0.0)).dot(matrix_w)
        return Y
    #
    # get matrix from data ...
    raw_matrix = dat.reset_index(drop=True)
    print raw_matrix
    # normalize the matrix ...
    normal_matrix = (raw_matrix - raw_matrix.mean())/raw_matrix.std()
    # get covariation matrix ...
    cov_matrix = normal_matrix.cov()
    # get matrix's eigen-vectors and values ...
    # BEWARE: use np.linalg.eigh, that assumes the symmetry of the matrix, to avoid imaginary numbers.
    eig_vals, eig_vecs = np.linalg.eigh(cov_matrix)
    print eig_vals
    # get varioation explained and sort everything by it ...
    fracs_var_explained, sorted_eigvecs = calculateVarianceExplainedSort(eig_vals, eig_vecs)
    # get projection matrix ...
    if matrix_w is None:
        matrix_w = getProjectionMatrix(sorted_eigvecs)
    else:
        # use matrix_w that was passed to the fiunction ...
        pass
    # Project data to the axes of highest variation (eig vectors)
    # dot product of PCA table and eigvecs ...
    Y = getDotProduct(normal_matrix, matrix_w)
    # in 'Y' - columns are Principal Components and rows correpond to sample ...
    # Y.shape -(num_rows,num_cols) ...
    Y_num_rows, Y_num_cols = Y.shape
    PC_dict = dict( ('PC%d'%(idx+1), Y[:,idx]) for idx in range(Y_num_cols) )
    var_dict = dict( ('PC%d'%(idx+1), frac) for idx,frac in enumerate(fracs_var_explained) )    
    return (var_dict,PC_dict,matrix_w)











def libMSA_extraction(reference_lib_fname,aln_fname,dataset=['Ss','Tm','Tt']):
    # get aln and indexes of templtes ...
    aln, tmp_idx_in_aln = readAlnIdentifyTemplates( aln_fname, template_ids=dataset )
    # extract positions we care about from 'reference_lib_fname' ...
    lib_dat = pd.read_csv(reference_lib_fname)
    # example  of '' organism-pos content: 'Ss-55' ...
    lib_dat['organism']  = lib_dat['organism-pos'].str.split('-').apply(lambda x: x[0])
    lib_dat['pos']       = lib_dat['organism-pos'].str.split('-').apply(lambda x: int(x[1]))
    # wt amino acid column MUST be in lib_dat!
    assert 'wtaa' in lib_dat.columns
    # separate lib_dat by organism ...
    tmp_dat_grouped = lib_dat.groupby('organism')
    # extract alignments corresponding to dataset ...
    #
    # empty storage ...
    lib_pos_dict    = {}
    matrix_list     = []
    aln_info_dict   = {}
    PC_dict_of_dict = {}
    for tmpid in dataset:
        tmp_index_aln = tmp_idx_in_aln[tmpid]
        tmplib_pos = tmp_dat_grouped.get_group(tmpid)['pos']
        lib_pos_dict[tmpid], lib_aas_in_aln, aln_for_lib = getAlnForTemplateLib(aln,tmp_index_aln,tmplib_pos,save_aln=True,sub_aln_fname='NOT_so_PRALINE499_%s.fasta'%tmpid)
        # do some echecking:
        # check that tmp sequence extracted from aln coincides with the library sequence ...
        print 'The following sequences must match, unless there is any trimming:'
        print ''.join(tmp_dat_grouped.get_group(tmpid)['wtaa']) 
        print ''.join(lib_aas_in_aln)
    #################################
    return



def MSA_pipeline(reference_lib_fname,aln_fname,dataset=['Ss','Tm','Tt'],verbose_out=False):
    # get aln and indexes of templtes ...
    aln, tmp_idx_in_aln = readAlnIdentifyTemplates( aln_fname, template_ids=dataset )
    # extract positions we care about from 'reference_lib_fname' ...
    lib_dat = pd.read_csv(reference_lib_fname)
    # example  of '' organism-pos content: 'Ss-55' ...
    lib_dat['organism']  = lib_dat['organism-pos'].str.split('-').apply(lambda x: x[0])
    lib_dat['pos']       = lib_dat['organism-pos'].str.split('-').apply(lambda x: int(x[1]))
    # wt amino acid column MUST be in lib_dat!
    assert 'wtaa' in lib_dat.columns
    # separate lib_dat by organism ...
    tmp_dat_grouped = lib_dat.groupby('organism')
    # extract alignments corresponding to dataset ...
    #
    # empty storage ...
    lib_pos_dict    = {}
    matrix_dict     = {}
    aln_info_dict   = {}
    matrix_w_dict   = {}
    PC_dict_of_dict = {}
    for tmpid in dataset:
        tmp_index_aln = tmp_idx_in_aln[tmpid]
        tmplib_pos = tmp_dat_grouped.get_group(tmpid)['pos']
        lib_pos_dict[tmpid], lib_aas_in_aln, aln_for_lib = getAlnForTemplateLib(aln,tmp_index_aln,tmplib_pos,save_aln=False,sub_aln_fname='')
        # do some echecking:
        # check that tmp sequence extracted from aln coincides with the library sequence ...
        print 'The following sequences must match, unless there is any trimming:'
        print ''.join(tmp_dat_grouped.get_group(tmpid)['wtaa']) 
        print ''.join(lib_aas_in_aln)
        #
        # get summary for dataset constituents independently ...
        aln_info_dict[tmpid] = getSummary(aln_for_lib)
        # get AA usage matrices ...
        # matrix_list.append( alnToAAmat(aln_for_lib, method='counts',pos_labels=tmp_dat_grouped.get_group(tmpid)['organism-pos'],log_transform=True) )
        matrix = alnToAAmat(aln_for_lib, method='counts',pos_labels=tmp_dat_grouped.get_group(tmpid)['organism-pos'],log_transform=True)
        # total_matrix = pd.concat( matrix_list )
        matrix_dict[tmpid] = matrix
        #
        # print matrix
        var_dict, PC_dict, matrix_w = runPCA(matrix)
        # print fraction of variability the PCs explain ...
        print
        print "variation explained by components for dataset: ",tmpid
        for pc in sorted( var_dict, key=lambda x: int(x.strip('PC')) ):
            print pc,'%.1f%%'%var_dict[pc]
        # store PCs ...
        PC_dict_of_dict[tmpid] = PC_dict
        matrix_w_dict[tmpid]   = matrix_w
    #
    lib_dat['wtKD'] = lib_dat['wtaa'].map(KD)
    #
    # turn aln summary info into a nicely labeled piece of DataFrame ...
    lib_dat = lib_dat.merge(pd.concat( pd.DataFrame(aln_info_dict[tmpid], index=tmp_dat_grouped.get_group(tmpid)['organism-pos']) for tmpid in dataset ),
        left_on='organism-pos',right_index=True)
    # now merge components to lib_dat as well ...
    # merge PC table to lib_dat as well ...
    lib_dat = lib_dat.merge( pd.concat( pd.DataFrame(PC_dict_of_dict[tmpid], index=tmp_dat_grouped.get_group(tmpid)['organism-pos']) for tmpid in dataset ),
        left_on='organism-pos', right_index=True )
    # return our output the large table ...
    if not verbose_out:
        return lib_dat, matrix_w_dict
    else:
        return lib_dat, matrix_w_dict, matrix_dict




def MSSA_pipeline(reference_lib_fname,aln_fname,dataset=['Ss','Tm','Tt'],verbose_out=False):
    # extract positions we care about from 'reference_lib_fname' ...
    lib_dat = pd.read_csv(reference_lib_fname)
    # example  of '' organism-pos content: 'Ss-55' ...
    lib_dat['organism']  = lib_dat['organism-pos'].str.split('-').apply(lambda x: x[0])
    lib_dat['pos']       = lib_dat['organism-pos'].str.split('-').apply(lambda x: int(x[1]))
    # wt amino acid column MUST be in lib_dat!
    assert 'wtaa' in lib_dat.columns
    # aln_fname must be a function in this case, just because of the way we've done structural alignments ...
    assert isfunction(aln_fname)
    # separate lib_dat by organism ...
    tmp_dat_grouped = lib_dat.groupby('organism')
    # extract alignments corresponding to dataset ...
    #
    # empty storage ...
    matrix_dict     = {}
    aln_info_dict   = {}
    matrix_w_dict   = {}
    PC_dict_of_dict = {}
    for tmpid in dataset:
        tmplib_pos = tmp_dat_grouped.get_group(tmpid)['pos']
        # get aln of templtes ...
        aln_for_lib, var_to_be_neglected = readAlnIdentifyTemplates(aln_fname(tmpid), template_ids=[] )
        #
        print "Checking is already done for structural alignments ..."
        # get summary for dataset constituents independently ...
        aln_info_dict[tmpid] = getSummary(aln_for_lib)
        # get AA usage matrices ...
        # matrix_list.append( alnToAAmat(aln_for_lib, method='counts',pos_labels=tmp_dat_grouped.get_group(tmpid)['organism-pos'],log_transform=True) )
        matrix = alnToAAmat(aln_for_lib, method='counts',pos_labels=tmp_dat_grouped.get_group(tmpid)['organism-pos'],log_transform=True)
        # total_matrix = pd.concat( matrix_list )
        matrix_dict[tmpid] = matrix
        #
        # print matrix
        var_dict, PC_dict, matrix_w = runPCA(matrix)
        # print fraction of variability the PCs explain ...
        print
        print "variation explained by components for dataset: ",tmpid
        for pc in sorted( var_dict, key=lambda x: int(x.strip('PC')) ):
            print pc,'%.1f%%'%var_dict[pc]
        # store PCs ...
        PC_dict_of_dict[tmpid] = PC_dict
        matrix_w_dict[tmpid]   = matrix_w
    #
    lib_dat['wtKD'] = lib_dat['wtaa'].map(KD)
    #
    # turn aln summary info into a nicely labeled piece of DataFrame ...
    lib_dat = lib_dat.merge(pd.concat( pd.DataFrame(aln_info_dict[tmpid], index=tmp_dat_grouped.get_group(tmpid)['organism-pos']) for tmpid in dataset ),
        left_on='organism-pos',right_index=True)
    # now merge components to lib_dat as well ...
    # merge PC table to lib_dat as well ...
    lib_dat = lib_dat.merge( pd.concat( pd.DataFrame(PC_dict_of_dict[tmpid], index=tmp_dat_grouped.get_group(tmpid)['organism-pos']) for tmpid in dataset ),
        left_on='organism-pos', right_index=True )
    # return our output the large table ...
    if not verbose_out:
        return lib_dat, matrix_w_dict
    else:
        return lib_dat, matrix_w_dict, matrix_dict





def MSA_pipeline_INTERSECT(reference_lib_fname,aln_fname,dataset=['Ss','Tm','Tt']):
    # get aln and indexes of templtes ...
    aln, tmp_idx_in_aln = readAlnIdentifyTemplates( aln_fname, template_ids=dataset )
    # extract positions we care about from 'reference_lib_fname' ...
    lib_dat = pd.read_csv(reference_lib_fname)
    # example  of '' organism-pos content: 'Ss-55' ...
    lib_dat['organism']  = lib_dat['organism-pos'].str.split('-').apply(lambda x: x[0])
    lib_dat['pos']       = lib_dat['organism-pos'].str.split('-').apply(lambda x: int(x[1]))
    # wt amino acid column MUST be in lib_dat!
    assert 'wtaa' in lib_dat.columns
    # separate lib_dat by organism ...
    tmp_dat_grouped = lib_dat.groupby('organism')
    # extract alignments corresponding to dataset ...
    #
    lib_pos_dict    = {}
    matrix_list     = []
    aln_info_dict   = {}
    matrix_w_dict   = {}
    PC_dict_of_dict = {}
    intersection_alignments = {}
    for tmpid in dataset:
        tmp_index_aln = tmp_idx_in_aln[tmpid]
        tmplib_pos = tmp_dat_grouped.get_group(tmpid)['pos']
        lib_pos_dict[tmpid], lib_aas_in_aln, aln_for_lib = getAlnForTemplateLib(aln,tmp_index_aln,tmplib_pos,save_aln=False,sub_aln_fname='')
        intersection_alignments[tmpid] = aln_for_lib
        # do some echecking:    # empty storage ...
    # FIND INTERSECTION BETWEEN ALL lib_pos_dict lists ...
    lib_pos_df = pd.DataFrame(lib_pos_dict)
    lib_pos_df['match'] = lib_pos_df.apply( lambda row: row.unique().size==1, axis=1 )
    pos_intersect_index = lib_pos_df[ lib_pos_df['match'] ].index
    print
    print "The size of interestion alignment is %d"%pos_intersect_index.size
    print
    #
    ################################
    # get sub-alignment from aln, given a set of positions of interest ...
    def get_subalignment(positions,aln):
        positions = list(positions)
        pos1 = positions[0]
        sub_aln = aln[:,pos1:pos1+1]
        for pos in positions[1:]:
            sub_aln += aln[:,pos:pos+1]
        return sub_aln
    ################################
    # 
    ######################################
    # DO THIS ALL OVER AGAIN BUT WITH THE INTERSECTION ALIGNMENT ... 
    ######################################
    # this can probably be done just a single time!!!!!!!!!!!!!!!!!!!!
    # to be continued ....
    for tmpid in dataset:
        aln_for_lib = get_subalignment(pos_intersect_index,intersection_alignments[tmpid])
        # get summary for dataset constituents independently ...
        aln_info_dict[tmpid] = getSummary(aln_for_lib)
        # get AA usage matrices ...
        # matrix_list.append( alnToAAmat(aln_for_lib, method='counts',pos_labels=tmp_dat_grouped.get_group(tmpid)['organism-pos'],log_transform=True) )
        pos_labels = tmp_dat_grouped.get_group(tmpid).iloc[pos_intersect_index]['organism-pos'] 
        matrix = alnToAAmat(aln_for_lib, method='counts',pos_labels=pos_labels,log_transform=True)
        # total_matrix = pd.concat( matrix_list )
        #
        # print matrix
        var_dict, PC_dict, matrix_w = runPCA(matrix)
        # print fraction of variability the PCs explain ...
        print
        print "variation explained by components for dataset: ",tmpid
        for pc in sorted( var_dict, key=lambda x: int(x.strip('PC')) ):
            print pc,'%.1f%%'%var_dict[pc]
        # store PCs ...
        PC_dict_of_dict[tmpid] = PC_dict
        matrix_w_dict[tmpid]   = matrix_w
    #
    lib_dat['wtKD'] = lib_dat['wtaa'].map(KD)
    #
    #
    #
    for tmpid in dataset:
        print pd.DataFrame(aln_info_dict[tmpid], index=tmp_dat_grouped.get_group(tmpid).iloc[pos_intersect_index]['organism-pos'])
    #
    for tmpid in dataset:
        print pd.DataFrame(PC_dict_of_dict[tmpid], index=tmp_dat_grouped.get_group(tmpid).iloc[pos_intersect_index]['organism-pos'])
    #
    # turn aln summary info into a nicely labeled piece of DataFrame ...
    lib_dat = lib_dat.merge(pd.concat( pd.DataFrame(aln_info_dict[tmpid], index=tmp_dat_grouped.get_group(tmpid).iloc[pos_intersect_index]['organism-pos']) for tmpid in dataset ),
        left_on='organism-pos',right_index=True,how='left')
    # now merge components to lib_dat as well ...
    # merge PC table to lib_dat as well ...
    lib_dat = lib_dat.merge( pd.concat( pd.DataFrame(PC_dict_of_dict[tmpid], index=tmp_dat_grouped.get_group(tmpid).iloc[pos_intersect_index]['organism-pos']) for tmpid in dataset ),
        left_on='organism-pos', right_index=True,how='left')
    # return our output the large table ...
    return lib_dat, matrix_w_dict



def MSA_pipeline_W(reference_lib_fname,proj_mat,aln_fname,dataset=['Ss','Tm','Tt']):
    # get aln and indexes of templtes ...
    aln, tmp_idx_in_aln = readAlnIdentifyTemplates( aln_fname, template_ids=dataset )
    # extract positions we care about from 'reference_lib_fname' ...
    lib_dat = pd.read_csv(reference_lib_fname)
    # example  of '' organism-pos content: 'Ss-55' ...
    lib_dat['organism']  = lib_dat['organism-pos'].str.split('-').apply(lambda x: x[0])
    lib_dat['pos']       = lib_dat['organism-pos'].str.split('-').apply(lambda x: int(x[1]))
    # wt amino acid column MUST be in lib_dat!
    assert 'wtaa' in lib_dat.columns
    # separate lib_dat by organism ...
    tmp_dat_grouped = lib_dat.groupby('organism')
    # extract alignments corresponding to dataset ...
    #
    # empty storage ...
    lib_pos_dict    = {}
    matrix_list     = []
    aln_info_dict   = {}
    PC_dict_of_dict = {}
    for tmpid in dataset:
        tmp_index_aln = tmp_idx_in_aln[tmpid]
        tmplib_pos = tmp_dat_grouped.get_group(tmpid)['pos']
        lib_pos_dict[tmpid], lib_aas_in_aln, aln_for_lib = getAlnForTemplateLib(aln,tmp_index_aln,tmplib_pos,save_aln=False,sub_aln_fname='')
        # do some echecking:
        # check that tmp sequence extracted from aln coincides with the library sequence ...
        print 'The following sequences must match, unless there is any trimming:'
        print ''.join(tmp_dat_grouped.get_group(tmpid)['wtaa']) 
        print ''.join(lib_aas_in_aln)
        #
        # get summary for dataset constituents independently ...
        aln_info_dict[tmpid] = getSummary(aln_for_lib)
        # get AA usage matrices ...
        # matrix_list.append( alnToAAmat(aln_for_lib, method='counts',pos_labels=tmp_dat_grouped.get_group(tmpid)['organism-pos'],log_transform=True) )
        matrix = alnToAAmat(aln_for_lib, method='counts',pos_labels=tmp_dat_grouped.get_group(tmpid)['organism-pos'],log_transform=True)
        # total_matrix = pd.concat( matrix_list )
        #
        print matrix
        var_dict, PC_dict, _ = runPCA(matrix,matrix_w=proj_mat)
        # print fraction of variability the PCs explain ...
        print
        print "variation explained by components for dataset: ",tmpid
        for pc in sorted( var_dict, key=lambda x: int(x.strip('PC')) ):
            print pc,'%.1f%%'%var_dict[pc]
        # store PCs ...
        PC_dict_of_dict[tmpid] = PC_dict
    #
    lib_dat['wtKD'] = lib_dat['wtaa'].map(KD)
    #
    # turn aln summary info into a nicely labeled piece of DataFrame ...
    lib_dat = lib_dat.merge(pd.concat( pd.DataFrame(aln_info_dict[tmpid], index=tmp_dat_grouped.get_group(tmpid)['organism-pos']) for tmpid in dataset ),
        left_on='organism-pos',right_index=True)
    # now merge components to lib_dat as well ...
    # merge PC table to lib_dat as well ...
    lib_dat = lib_dat.merge( pd.concat( pd.DataFrame(PC_dict_of_dict[tmpid], index=tmp_dat_grouped.get_group(tmpid)['organism-pos']) for tmpid in dataset ),
        left_on='organism-pos', right_index=True )
    # return our output the large table ...
    return lib_dat




def MSSA_pipeline_W(reference_lib_fname,proj_mat,aln_fname,dataset=['Ss','Tm','Tt']):
    # extract positions we care about from 'reference_lib_fname' ...
    lib_dat = pd.read_csv(reference_lib_fname)
    # example  of '' organism-pos content: 'Ss-55' ...
    lib_dat['organism']  = lib_dat['organism-pos'].str.split('-').apply(lambda x: x[0])
    lib_dat['pos']       = lib_dat['organism-pos'].str.split('-').apply(lambda x: int(x[1]))
    # wt amino acid column MUST be in lib_dat!
    assert 'wtaa' in lib_dat.columns
    # separate lib_dat by organism ...
    tmp_dat_grouped = lib_dat.groupby('organism')
    # extract alignments corresponding to dataset ...
    #
    # empty storage ...
    matrix_list     = []
    aln_info_dict   = {}
    PC_dict_of_dict = {}
    for tmpid in dataset:
        tmplib_pos = tmp_dat_grouped.get_group(tmpid)['pos']
        # get aln of templtes ...
        aln_for_lib, var_to_be_neglected = readAlnIdentifyTemplates('../SPalignNS/%s_tim_aln.fasta'%tmpid, template_ids=[] )
        #
        print "Checking is already done for structural alignments ..."
        # get summary for dataset constituents independently ...
        aln_info_dict[tmpid] = getSummary(aln_for_lib)
        # get AA usage matrices ...
        # matrix_list.append( alnToAAmat(aln_for_lib, method='counts',pos_labels=tmp_dat_grouped.get_group(tmpid)['organism-pos'],log_transform=True) )
        matrix = alnToAAmat(aln_for_lib, method='counts',pos_labels=tmp_dat_grouped.get_group(tmpid)['organism-pos'],log_transform=True)
        # total_matrix = pd.concat( matrix_list )
        #
        # print matrix
        var_dict, PC_dict, _ = runPCA(matrix,matrix_w=proj_mat)
        # print fraction of variability the PCs explain ...
        print
        print "variation explained by components for dataset: ",tmpid
        for pc in sorted( var_dict, key=lambda x: int(x.strip('PC')) ):
            print pc,'%.1f%%'%var_dict[pc]
        # store PCs ...
        PC_dict_of_dict[tmpid] = PC_dict
    #
    lib_dat['wtKD'] = lib_dat['wtaa'].map(KD)
    #
    # turn aln summary info into a nicely labeled piece of DataFrame ...
    lib_dat = lib_dat.merge(pd.concat( pd.DataFrame(aln_info_dict[tmpid], index=tmp_dat_grouped.get_group(tmpid)['organism-pos']) for tmpid in dataset ),
        left_on='organism-pos',right_index=True)
    # now merge components to lib_dat as well ...
    # merge PC table to lib_dat as well ...
    lib_dat = lib_dat.merge( pd.concat( pd.DataFrame(PC_dict_of_dict[tmpid], index=tmp_dat_grouped.get_group(tmpid)['organism-pos']) for tmpid in dataset ),
        left_on='organism-pos', right_index=True )
    # return our output the large table ...
    return lib_dat









# http://stackoverflow.com/questions/419163/what-does-if-name-main-do
if __name__ == "__main__":
    print "main program body ..."
    #
    #
    #
    feat_fname = '../AccessDB-files/features-original.csv'
    # data_fname = '../AccessDB-files/db-fitness.csv'
    # msa_fname  = '../AccessDB-files/IGPS_seq_Praline_align/new-PRALINE-aln.fasta'
    # msa_fname  = './PFam/IGPS-Pfam-Aln.txt'
    msa_fname  = './PFam/PF00218_full-Aln.txt'
    #
    #
    # libMSA_extraction(feat_fname,msa_fname,dataset=['Ss','Tm','Tt'])
    libMSA_extraction(feat_fname,msa_fname,dataset=['Ss','Tm'])
    #
    # pass
    # ######################################################################
    # #  MSA PRALINE output ...
    # ######################################################################
    # # Ss_dat['gap_pos']==Tm_dat['gap_pos']
    # # positions_chosen = optimist_pos
    # positions_chosen = conserv_pos
    # ################################################################
    # # get PC1 aligned to  the library ...
    # aligned_PC1 = copy.deepcopy(Ss_dat['gap_pos'])
    # aligned_PC1[:] = None
    # aligned_PC1.loc[positions_chosen.index] = Y[:,0]
    # print aligned_PC1.head()
    # aligned_PC1.to_csv('IGPS_MSA_PC1.csv', index=True, sep=',', na_rep='')
    # ################################################################
    # # get PC2 aligned to  the library ...
    # aligned_PC2 = copy.deepcopy(Ss_dat['gap_pos'])
    # aligned_PC2[:] = None
    # aligned_PC2.loc[positions_chosen.index] = Y[:,1]
    # print aligned_PC2.head()
    # aligned_PC2.to_csv('IGPS_MSA_PC2.csv', index=True, sep=',', na_rep='')
    # ################################################################
    # # get PC2 aligned to  the library ...
    # aligned_ic = copy.deepcopy(Ss_dat['gap_pos'])
    # aligned_ic[:] = None
    # aligned_ic.loc[positions_chosen.index] = ic_vector
    # print aligned_ic.head()
    # aligned_ic.to_csv('IGPS_ic_vector.csv', index=True, sep=',', na_rep='')















































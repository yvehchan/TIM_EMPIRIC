
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats as st
import matplotlib.mlab as mlab
from Bio import SeqUtils
from Bio.SeqUtils import ProtParamData 
import math
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import collections


# EMPIRIC_raw_data_fname = "db-fitness.csv"
# # most important reference file with all the libraries information ...
# EMPIRIC_features_fname = "features-original.csv"
# 20 amino acids ...
aacids = sorted(list(SeqUtils.IUPAC.protein.letters))


#filter dataset 
def filterDataset(dat, dataset):
    """This function filters low-qual and other controls data from the
    original expanded version of the EMPIRIC dataset."""
    #
    dat = dat[dat['organism'].isin(dataset)]
    no_mmei_index = dat['mmei']=='no'
    nonstop_index = dat['mutstop']=='no'
    zerofit_index = dat['fitness'].abs()>1e-4
    mutwt_index = dat['mutwt']=='no'
    dat = dat[no_mmei_index & nonstop_index & zerofit_index & mutwt_index]
    #print "Filtered data"
    return dat


#add features to table
def pivotAddFeatures(dat, dat_features):
    """Takes raw EMPIRIC data (unrolled and filtered), and features of positions in library:
    Returns a merged table, where each library position is characterized by features and 20-EMPIRIC fitness readings.
    BEWARE:
    Unique labels are 'organism-pos' and their order as in features table, which is cruicial for compatibility."""
    #pivot EMPIRIC dataset to have mutaa(fitness readings) as columns and wt pos as rows
    table = dat.pivot(index='organism-pos',columns='mutaa',values='fitness')
    # add wtaa residues to table based on org-pos
    columns_to_merge = dat[['organism-pos','wtaa']].drop_duplicates()
    merged_table = table.merge(columns_to_merge, how='inner',left_index=True,right_on='organism-pos').reset_index(drop=True)
    # add other features of interest to the table
    # MUST merge the table to features to keep org-pos order as in EMPIRIC_features_fname ...
    merged_table = dat_features.merge(merged_table,how='inner',on='organism-pos',suffixes=('', '_copy')).reset_index(drop=True)
    #add kd of wildtype
    KD = ProtParamData.kd
    merged_table['wtKD'] = merged_table['wtaa'].map(KD)
    # we can add more features later on here ...
    # ...
    # ...
    #print merged_table
    return merged_table





#main function for running PCA, calls on subfunctions
def runPCA(dat):
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
    # normalize the matrix ...
    normal_matrix = (raw_matrix - raw_matrix.mean())/raw_matrix.std()
    # get covariation matrix ...
    cov_matrix = normal_matrix.cov()
    # get matrix's eigen-vectors and values ...
    # BEWARE: use np.linalg.eigh, that assumes the symmetry of the matrix, to avoid imaginary numbers.
    eig_vals, eig_vecs = np.linalg.eigh(cov_matrix)
    # get varioation explained and sort everything by it ...
    fracs_var_explained, sorted_eigvecs = calculateVarianceExplainedSort(eig_vals, eig_vecs)
    # get projection matrix ...
    matrix_w = getProjectionMatrix(sorted_eigvecs)  
    # Project data to the axes of highest variation (eig vectors)
    # dot product of PCA table and eigvecs ...
    Y = getDotProduct(normal_matrix, matrix_w)
    # in 'Y' - columns are Principal Components and rows correpond to sample ...
    # Y.shape -(num_rows,num_cols) ...
    Y_num_rows, Y_num_cols = Y.shape
    PC_dict = dict( ('PC%d'%(idx+1), Y[:,idx]) for idx in range(Y_num_cols) )
    var_dict = dict( ('PC%d'%(idx+1), frac) for idx,frac in enumerate(fracs_var_explained) )    
    return (var_dict,PC_dict,matrix_w)







def EMPIRIC_pipeline_separate(EMPIRIC_features_fname, EMPIRIC_raw_data_fname, dataset=['Ss','Tm','Tt']):
    # read raw data ...
    dat = pd.read_csv(EMPIRIC_raw_data_fname)
    # filter and extract dataset of interest ...
    dat = filterDataset(dat, dataset=dataset)
    # extract positions we care about from 'EMPIRIC_features_fname' ...
    lib_dat = pd.read_csv(EMPIRIC_features_fname)
    # example  of '' organism-pos content: 'Ss-55' ...
    lib_dat['organism']  = lib_dat['organism-pos'].str.split('-').apply(lambda x: x[0])
    lib_dat['pos']       = lib_dat['organism-pos'].str.split('-').apply(lambda x: x[1])
    # wt amino acid column MUST be in lib_dat!
    assert 'wtaa' in lib_dat.columns
    #
    # merge pivoted raw data with the features ...
    merged_dat = pivotAddFeatures(dat, lib_dat)
    #
    # separate lib_dat by organism ...
    tmp_dat_grouped = lib_dat.groupby('organism')
    # extract alignments corresponding to dataset ...
    # empty storage ...
    lib_pos_dict    = {}
    matrix_list     = []
    matrix_w_dict   = {}
    aln_info_dict   = {}
    PC_dict_of_dict = {}
    for tmpid in dataset:
        print
        print "Preparing separate PCAs for different organisms: %s"%tmpid
        matrix = merged_dat[merged_dat['organism']==tmpid][aacids]
        # print matrix
        var_dict, PC_dict, matrix_w = runPCA(matrix)
        # print fraction of variability the PCs explain ...
        print
        print "variation explained by components for dataset: ",tmpid
        for pc in sorted( var_dict, key=lambda x: int(x.strip('PC')) ):
            print pc,'%.1f%%'%var_dict[pc]
        # store PCs and projection matrices ...
        PC_dict_of_dict[tmpid] = PC_dict
        matrix_w_dict[tmpid]   = matrix_w
    # # PCA ... 
    #
    # now merge components to lib_dat as well ...
    # merge PC table to lib_dat as well ...
    merged_dat = merged_dat.merge(
                    pd.concat(pd.DataFrame( PC_dict_of_dict[tmpid],index=merged_dat[merged_dat['organism']==tmpid]['organism-pos']) for tmpid in dataset),
                    left_on='organism-pos', right_index=True )
    # add more columns to the merged_dat ...
    # get average fitness per position, along with min/max fitness ...
    merged_dat['fitness']     = merged_dat[aacids].mean(axis=1)
    merged_dat['min_fitness'] = merged_dat[aacids].min(axis=1)
    merged_dat['max_fitness'] = merged_dat[aacids].max(axis=1)
    # return our output the large table ...
    return merged_dat, matrix_w_dict






def EMPIRIC_pipeline_joint(EMPIRIC_features_fname, EMPIRIC_raw_data_fname, dataset=['Ss','Tm','Tt']):
    # read raw data ...
    dat = pd.read_csv(EMPIRIC_raw_data_fname)
    # filter and extract dataset of interest ...
    dat = filterDataset(dat, dataset=dataset)
    # extract positions we care about from 'EMPIRIC_features_fname' ...
    lib_dat = pd.read_csv(EMPIRIC_features_fname)
    # example  of '' organism-pos content: 'Ss-55' ...
    lib_dat['organism']  = lib_dat['organism-pos'].str.split('-').apply(lambda x: x[0])
    lib_dat['pos']       = lib_dat['organism-pos'].str.split('-').apply(lambda x: x[1])
    # wt amino acid column MUST be in lib_dat!
    assert 'wtaa' in lib_dat.columns
    #
    # merge pivoted raw data with the features ...
    merged_dat = pivotAddFeatures(dat, lib_dat)
    #
    # separate lib_dat by organism ...
    tmp_dat_grouped = lib_dat.groupby('organism')
    # extract alignments corresponding to dataset ...
    # PCA ... 
    total_matrix = merged_dat[aacids] 
    var_dict, PC_dict, matrix_w = runPCA(total_matrix)
    # print fraction of variability the PCs explain ...
    print
    print "variation explained by components for dataset: ",dataset
    for pc in sorted( var_dict, key=lambda x: int(x.strip('PC')) ):
        print pc,'%.1f%%'%var_dict[pc]
    #
    # now merge components to lib_dat as well ...
    # merge PC table to lib_dat as well ...
    merged_dat = merged_dat.merge( pd.DataFrame(PC_dict,index=merged_dat['organism-pos']), left_on='organism-pos', right_index=True )
    # add more columns to the merged_dat ...
    # get average fitness per position, along with min/max fitness ...
    merged_dat['fitness']     = merged_dat[aacids].mean(axis=1)
    merged_dat['min_fitness'] = merged_dat[aacids].min(axis=1)
    merged_dat['max_fitness'] = merged_dat[aacids].max(axis=1)
    # return our output the large table ...
    return merged_dat, matrix_w













# http://stackoverflow.com/questions/419163/what-does-if-name-main-do
if __name__ == "__main__":
    print "main program body ..."
    pass














































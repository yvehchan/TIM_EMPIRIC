
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


EMPIRIC_raw_data_fname = "db-fitness.csv"
# most important reference file with all the libraries information ...
EMPIRIC_features_fname = "features-original.csv"
# 20 amino acids ...
aacids = list(SeqUtils.IUPAC.protein.letters)


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
    #print dat
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
    merged_table = dat_features.merge(merged_table,how='inner',on='organism-pos').reset_index(drop=True)
    #add kd of wildtype
    KD = ProtParamData.kd
    merged_table['wtKD'] = merged_table['wtaa'].map(KD)
    # we can add more features later on here ...
    # ...
    # ...
    #print merged_table
    return merged_table



#main function for running PCA, calls on subfunctions
def runPCA(dat, num_eigenv, pc1, pc2, filtParam):
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
        fracs_var_explained = [eigval/total_eigv*100.0 for eigval,eigvec in sorted_eigs]
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
    matrix_w = getMatrix(sorted_eigs, num_eigenv)  
    # Project data to the axes of highest variation (eig vectors)
    # dot product of PCA table and eigvecs ...
    Y = getDotProduct(normal_matrix, matrix_w)
    # in 'Y' - columns are Principal Components and rows correpond to sample ...
    # Y.shape -(num_rows,num_cols) ...
    Y_num_rows, Y_num_cols = Y.shape
    PC_dict = dict( ('PC%d'%(idx+1), Y[:,idx]) for idx in range(Y_num_cols) )
    var_dict = dict( ('PC%d'%(idx+1), frac) for idx,frac in enumerate(fracs_var_explained) )    
    return (var_dict,PC_dict)




# http://stackoverflow.com/questions/419163/what-does-if-name-main-do
if __name__ == "__main__":
    print "main program body ..."
    pass














































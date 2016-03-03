
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










































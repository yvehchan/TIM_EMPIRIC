import EMPIRIC_analysis_EXP as emp
import MSA_analysis_EXP as msa
import pandas  as pd
from scipy import stats as st
# !pwd
import matplotlib.pyplot as plt
#!ls
import numpy as np
from Bio import SeqUtils
####################
import matplotlib.mlab as mlab
from Bio.SeqUtils import ProtParamData 
import math
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import collections
import copy
#import seaborn as sns

feat_fname = './features-original.csv'
data_fname = './db-fitness.csv'

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
	##allow stops and wt
    dat = dat[no_mmei_index]
#    dat = dat[no_mmei_index & nonstop_index & zerofit_index & mutwt_index]
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

    return merged_table



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
    return merged_dat


dat = EMPIRIC_pipeline_joint(feat_fname, data_fname, dataset=['Ss','Tm','Tt'])



###################################
#  Where is the border between benifitial/delitirious, do one need the wt-like as a separate???
###################################

# (1)
# solution one (simplest) chose fit=0.0, as separation between 2 groups: deleterious AND benifitial
fit_threshold = 0.0


# (2)
# choose threshold as the fitness value that corresponds to the local minimum of the saddle on the overall bimodal fitness distribution.
flat_fit_dat = dat[aacids].unstack()
flat_fit_dat = flat_fit_dat.reset_index(drop=True)
flat_fit_dat = flat_fit_dat[flat_fit_dat.notnull()]
# get histogramm ...
# ACHTUNG!!! value to be determined is unstable with restect to the choice of the bins number.
h,b = np.histogram(flat_fit_dat,bins=300)
fit_hist = pd.DataFrame([h,b]).T
fit_hist = fit_hist.rename(columns={0:'count',1:'fit_bins'})
# consider saddle only ...
fit_hist = fit_hist[(fit_hist['fit_bins']>-1.0)&(fit_hist['fit_bins']<0.0)].reset_index(drop=True)
# finally choose the threshold fit value ...
fit_threshold = fit_hist['fit_bins'][fit_hist['count'].argmin()]

# (3)
# use approximartion of the bimodal distro with 2 gauss-like distros
# dunno how.

# (4)
# any ideas? common sense ranges for benifitial and deleterious ???


# function that would be getting ratio ...
def get_ratio(fff,threshold=fit_threshold):
    # fff has to be DataFrame with fitness values ...
    flat_fit_dat = fff[aacids].unstack().reset_index(drop=True)
    flat_fit_dat = flat_fit_dat[flat_fit_dat.notnull()]
    # get ratio finally ...
    num_benefit = (flat_fit_dat>threshold).sum()
    return num_benefit*100.0/flat_fit_dat.size



#####################################################################################################
# SETTING P-VALUES ON COMPARISON STATEMENTS ...
#####################################################################################################

#########################################################
#  Comparing 1 feature ...
########################################################

# null hypoth is that mean in different categories is the same
# from visual stand point implies that OBSERVED delta-mena lies "deep" inside the distribution.
def getp(feature, cat1,cat2):
  fff = copy.deepcopy(dat[aacids+[feature]])
  #filter by categories
  fff = fff[fff[feature].isin([cat1,cat2])]
  #getmean
  fff['l1'] = fff[feature]
  means_wt = fff.groupby('l1')[aacids].mean().mean(axis=1)
  #print means_wt
  #getrandomlist
  means = []
  for i in range(10000):
	fff['l1'] = np.random.permutation(fff[feature])
	means.append( fff.groupby('l1')[aacids].mean().mean(axis=1) )
  res = pd.DataFrame(map(pd.Series,means))
  #fff, res = getrandom(fff, feature, cat1, cat2)
  supporting_null = ( np.abs(res[cat1]-res[cat2]) >= np.abs(means_wt[cat1]-means_wt[cat2]) ).sum()
  total,_ = res.shape
  pval = float(supporting_null)/float(total)
  print "Comparing %s and %s yields Pval=%.5f\n"%(str(cat1),str(cat2),pval)
  return pval


# null hypoth is that fractrion of benefitial mutations in different categories is the same
# from visual stand point implies that OBSERVED delta-fraction lies "deep" inside the distribution.
def getp_ratio(feature, cat1,cat2):
  fff = copy.deepcopy(dat[aacids+[feature]])
  #filter by categories
  fff = fff[fff[feature].isin([cat1,cat2])]
  #getmean
  fff['l1'] = fff[feature]
  ratios_wt = fff.groupby('l1').apply(get_ratio)
  #print means_wt
  #getrandomlist
  ratios = []
  for i in range(10000):
    fff['l1'] = np.random.permutation(fff[feature])
    ratios.append( fff.groupby('l1').apply(get_ratio) )
  res = pd.DataFrame(map(pd.Series,ratios))
  #fff, res = getrandom(fff, feature, cat1, cat2)
  supporting_null = ( np.abs(res[cat1]-res[cat2]) >= np.abs(ratios_wt[cat1]-ratios_wt[cat2]) ).sum()
  total,_ = res.shape
  pval = float(supporting_null)/float(total)
  print "(RATIOS!) Comparing %s and %s yields Pval=%.5f\n"%(str(cat1),str(cat2),pval)
  return pval

  
#######################

#set fitness threshold
fit_threshold = 0.0


### 0=abloops, 1=baloops, 2=bstrands
print 
print "Secondary structure: "
getp("Loop-Type", 0,2)
getp("Loop-Type", 2,1)
getp("Loop-Type", 0,1)

getp_ratio("Loop-Type", 0,2)
getp_ratio("Loop-Type", 2,1)
getp_ratio("Loop-Type", 0,1)

#------bstrands in vs out 
print 
print "Comparing b-strands: In and out ..."
getp("sc_direction", "In", "Out")
getp_ratio("sc_direction", "In", "Out")


#########################################################
#  Comparing 2 features ...
########################################################
def getp(feature1, feature2, f1c1, f1c2, f2c1, f2c2):
  fff = copy.deepcopy(dat[aacids+[feature1,feature2]])
  #filter by categories
  fff = fff[fff[feature1].isin([f1c1,f1c2])]
  fff['l1'] = fff[feature1]
  fff = fff[fff[feature2].isin([f2c1,f2c2])]
  fff['l2'] = fff[feature2]
  #getmean
  means_wt = fff.groupby(['l1','l2'])[aacids].mean().mean(axis=1)
  #print means_wt
  means = []
  #getrandomlist
  for i in range(1000):
    a,b = np.random.permutation(fff[['l1','l2']]).T
    fff['l1'] = a
    fff['l2'] = b
    means.append( fff.groupby(['l1','l2'])[aacids].mean().mean(axis=1) )
  res = pd.DataFrame(map(pd.Series,means))
  supporting_null = ( np.abs(res[(f1c1,f2c1)]-res[(f1c2,f2c2)]) >= np.abs(means_wt[(f1c1,f2c1)]-means_wt[(f1c2,f2c2)]) ).sum()
  total,_ = res.shape
  pval = float(supporting_null)/float(total)
  print "Comparing  %s and  %s yields Pval=%.5f\n"%(str((f1c1,f2c1)),str((f1c2,f2c2)),pval)
  return pval


def getp_ratio(feature1, feature2, f1c1, f1c2, f2c1, f2c2):
  fff = copy.deepcopy(dat[aacids+[feature1,feature2]])
  #filter by categories
  fff = fff[fff[feature1].isin([f1c1,f1c2])]
  fff['l1'] = fff[feature1]
  fff = fff[fff[feature2].isin([f2c1,f2c2])]
  fff['l2'] = fff[feature2]
  #getmean
  ratios_wt = fff.groupby(['l1','l2']).apply(get_ratio)
  #print ratios_wt
  ratios = []
  #getrandomlist
  for i in range(1000):
    a,b = np.random.permutation(fff[['l1','l2']]).T
    fff['l1'] = a
    fff['l2'] = b
    ratios.append( fff.groupby(['l1','l2']).apply(get_ratio) )
  res = pd.DataFrame(map(pd.Series,ratios))
  supporting_null = ( np.abs(res[(f1c1,f2c1)]-res[(f1c2,f2c2)]) >= np.abs(ratios_wt[(f1c1,f2c1)]-ratios_wt[(f1c2,f2c2)]) ).sum()
  total,_ = res.shape
  pval = float(supporting_null)/float(total)
  print "(RATIOS!) Comparing  %s and  %s yields Pval=%.5f\n"%(str((f1c1,f2c1)),str((f1c2,f2c2)),pval)
  return pval



#example  
#feature1 = 'sslayers'
#feature2 = 'sc_direction'
#f1c1 = "Layer1"
#f1c2 = "Layer2"
#f2c1 = "In"
#f2c2 = "In"
#getp(feature1, feature2, f1c1, f1c2, f2c1, f2c2)

#######################
### categories 
print
print "Comparing ab-loop by position: even and odd"
getp("sc_direction", "StrandParity", "abloop", "abloop", "Odd", "Even")
getp_ratio("sc_direction", "StrandParity", "abloop", "abloop", "Odd", "Even")

print
print "Comparing ab-loop by position: even and odd"
getp("sslayers", "StrandParity", "ab-3", "ab-3", "Odd", "Even")
getp("sslayers", "StrandParity", "ab-2", "ab-2", "Odd", "Even")
getp("sslayers", "StrandParity", "ab-1", "ab-1", "Odd", "Even")

getp("sslayers", "StrandParity", "ab-3", "ab-2", "Odd", "Odd")
getp("sslayers", "StrandParity", "ab-3", "ab-1", "Odd", "Odd")
getp("sslayers", "StrandParity", "ab-2", "ab-1", "Odd", "Odd")

getp("sslayers", "StrandParity", "ab-3", "ab-2", "Even", "Even")
getp("sslayers", "StrandParity", "ab-3", "ab-1", "Even", "Even")
getp("sslayers", "StrandParity", "ab-2", "ab-1", "Even", "Even")

getp_ratio("sslayers", "StrandParity", "ab-3", "ab-3", "Odd", "Even")
getp_ratio("sslayers", "StrandParity", "ab-2", "ab-2", "Odd", "Even")
getp_ratio("sslayers", "StrandParity", "ab-1", "ab-1", "Odd", "Even")

getp_ratio("sslayers", "StrandParity", "ab-3", "ab-2", "Odd", "Odd")
getp_ratio("sslayers", "StrandParity", "ab-3", "ab-1", "Odd", "Odd")
getp_ratio("sslayers", "StrandParity", "ab-2", "ab-1", "Odd", "Odd")

getp_ratio("sslayers", "StrandParity", "ab-3", "ab-2", "Even", "Even")
getp_ratio("sslayers", "StrandParity", "ab-3", "ab-1", "Even", "Even")
getp_ratio("sslayers", "StrandParity", "ab-2", "ab-1", "Even", "Even")

print
print "Comparing bstrands by layer: in vs out"
getp("sc_direction", "sslayers", "In", "Out", "Layer1", "Layer1")
getp("sc_direction", "sslayers", "In", "Out", "Layer2", "Layer2")
getp("sc_direction", "sslayers", "In", "Out", "Layer3", "Layer3")
getp("sc_direction", "sslayers", "In", "Out", "Layer4", "Layer4")

getp_ratio("sc_direction", "sslayers", "In", "Out", "Layer1", "Layer1")
getp_ratio("sc_direction", "sslayers", "In", "Out", "Layer2", "Layer2")
getp_ratio("sc_direction", "sslayers", "In", "Out", "Layer3", "Layer3")
getp_ratio("sc_direction", "sslayers", "In", "Out", "Layer4", "Layer4")

print
print "Comparing bstrands (in): layers"
getp("sc_direction", "sslayers", "In", "In", "Layer1", "Layer2")
getp("sc_direction", "sslayers", "In", "In", "Layer1", "Layer3")
getp("sc_direction", "sslayers", "In", "In", "Layer1", "Layer4")
getp("sc_direction", "sslayers", "In", "In", "Layer2", "Layer3")
getp("sc_direction", "sslayers", "In", "In", "Layer2", "Layer4")
getp("sc_direction", "sslayers", "In", "In", "Layer3", "Layer4")

getp_ratio("sc_direction", "sslayers", "In", "In", "Layer1", "Layer2")
getp_ratio("sc_direction", "sslayers", "In", "In", "Layer1", "Layer3")
getp_ratio("sc_direction", "sslayers", "In", "In", "Layer1", "Layer4")
getp_ratio("sc_direction", "sslayers", "In", "In", "Layer2", "Layer3")
getp_ratio("sc_direction", "sslayers", "In", "In", "Layer2", "Layer4")
getp_ratio("sc_direction", "sslayers", "In", "In", "Layer3", "Layer4")


print
print "Comparing bstrands (out): layers"
getp("sc_direction", "sslayers", "Out", "Out", "Layer1", "Layer2")
getp("sc_direction", "sslayers", "Out", "Out", "Layer1", "Layer2")
getp("sc_direction", "sslayers", "Out", "Out", "Layer1", "Layer4")
getp("sc_direction", "sslayers", "Out", "Out", "Layer2", "Layer3")
getp("sc_direction", "sslayers", "Out", "Out", "Layer2", "Layer4")
getp("sc_direction", "sslayers", "Out", "Out", "Layer3", "Layer4")

getp_ratio("sc_direction", "sslayers", "Out", "Out", "Layer1", "Layer2")
getp_ratio("sc_direction", "sslayers", "Out", "Out", "Layer1", "Layer2")
getp_ratio("sc_direction", "sslayers", "Out", "Out", "Layer1", "Layer4")
getp_ratio("sc_direction", "sslayers", "Out", "Out", "Layer2", "Layer3")
getp_ratio("sc_direction", "sslayers", "Out", "Out", "Layer2", "Layer4")
getp_ratio("sc_direction", "sslayers", "Out", "Out", "Layer3", "Layer4")





















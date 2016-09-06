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

#import seaborn as sns

feat_fname = './Inputs/features-original.csv'
data_fname = './Inputs/db-fitness.csv'
version='includeWtandNull_biochem_res'
# 20 amino acids ...
aacids = sorted(list(SeqUtils.IUPAC.protein.letters))
#aabiochem= {'Basic','Acidic','Amidic','SulfurContaining','Hydroxylic','SmallAliphatic','Aromatic','BranchedAliphatic'}

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
cors_dict = {}
#prints out list of column headers
#list(dat.columns.values)

# ##########################################
# # corrs averaging by the same wt AA ...
# ###########################################

avg_corrs = []
all_corrs = []

for aa in aacids:
   yyy = dat[aacids][dat['wtaa']==aa]
   cmat = yyy.T.corr('pearson')
   n,n = cmat.shape
   tril_values = cmat.as_matrix()[np.tril_indices(n,-1)]
   avg_corr = tril_values.mean()
   avg_corrs.append((aa,avg_corr,n,0.5*n*(n-1)))
   # print aa,avg_corr
   all_corrs.append(tril_values)
   
all_corrs = np.hstack(all_corrs)
avg_corrs_df = pd.DataFrame(avg_corrs,columns=['aa','r','num','num_triangle'])
avg_corrs_df.to_csv('./numbers_for_paper/avg_corrs_ident_AA_%s.csv' %(version),index=False)
print "Median corr(R) over all identical AAs: ", np.median(all_corrs)
print "Mean corr(R) over all identical AAs: ", all_corrs.mean()
print "Mean corr(R) over all identical AAs(not weighted): ", avg_corrs_df['r'].mean()

with open('./numbers_for_paper/avg_corrs_ident_AA_%s.csv' %(version),'a') as fp:
   fp.write('\nmedian weighted, %.4f\nmean weighted, %.4f\nmean non-weighted, %.4f'%(np.median(all_corrs), all_corrs.mean(),avg_corrs_df['r'].mean()))

cors_dict['identical_aa'] = pd.Series(all_corrs)

 ##########################################
# # corrs averaging by the similar biochem AA ...
# ###########################################


avg_corrs = []
all_corrs = []
for feat in set(dat['biochem']):
 
   yyy = dat[aacids][dat['biochem']==feat]
   cmat = yyy.T.corr('pearson')
   n,n = cmat.shape
   tril_values = cmat.as_matrix()[np.tril_indices(n,-1)]
   avg_corr = tril_values.mean()
   avg_corrs.append((feat,avg_corr,n,0.5*n*(n-1)))
   print feat,avg_corr
   all_corrs.append(tril_values)
   
all_corrs = np.hstack(all_corrs)
avg_corrs_df = pd.DataFrame(avg_corrs,columns=['Biochem','r','num','num_triangle'])
avg_corrs_df.to_csv('./numbers_for_paper/avg_corrs_similar_biochem_AA_%s.csv' %(version),index=False)
print "Median corr(R) over all similar biochem AAs: ", np.median(all_corrs)
print "Mean corr(R) over all similar biochem AAs: ", all_corrs.mean()
print "Mean corr(R) over all similar biochem AAs(not weighted): ", avg_corrs_df['r'].mean()

with open('./numbers_for_paper/avg_corrs_similar_biochem_AA_%s.csv' %(version),'a') as fp:
   fp.write('\nmedian weighted, %.4f\nmean weighted, %.4f\nmean non-weighted, %.4f'%(np.median(all_corrs), all_corrs.mean(),avg_corrs_df['r'].mean()))

cors_dict['similar_aa'] = pd.Series(all_corrs)


# ##########################################
# # corrs averaging by structurally identical positions ...
# ###########################################

avg_corrs = []
all_corrs = []
for pos in range(1,97):
 yyy = dat[aacids][dat['align_orth_1-96']==pos]
 sigmasq = np.square(yyy.stack().std())
 avgfit = yyy.stack().mean()
 cmat = yyy.T.corr('pearson')
 n,n = cmat.shape
 tril_values = cmat.as_matrix()[np.tril_indices(n,-1)]
 avg_corr = tril_values.mean()
 avg_corrs.append((pos,avg_corr,n,0.5*n*(n-1), sigmasq, avgfit))
 #print aa,avg_corr
 all_corrs.append(tril_values)

all_corrs = np.hstack(all_corrs)

avg_corrs_df = pd.DataFrame(avg_corrs,columns=['pos','r','num','num_triangle', "sigmasq", "avgfit"])
avg_corrs_df.to_csv('./numbers_for_paper/avg_corrs_ident_pos_%s.csv' %(version),index=False)
# # # avg_corrs_df
print "Median corr(R) over all identical structural positions: ", np.median(all_corrs)
print "Mean corr(R) over all identical structural positions: ", all_corrs.mean()
print "Mean corr(R) over all identical structural positions(not weighted): ", avg_corrs_df['r'].mean()

with open('./numbers_for_paper/avg_corrs_ident_pos_%s.csv' %(version),'a') as fp:
 fp.write('\nmedian weighted, %.4f\nmean weighted, %.4f\nmean non-weighted, %.4f'%(np.median(all_corrs), all_corrs.mean(),avg_corrs_df['r'].mean()))


cors_dict['aligned_pos'] = pd.Series(all_corrs)



# ##########################################
# # corrs averaging by structurally identical positions AND 3 different WT residues...
# ###########################################

avg_corrs = []
all_corrs = []
wtlist = []

for pos in range(1,97):
 yyy = dat[aacids][dat['align_orth_1-96']==pos]
 wtlist = dat["wtaa"][dat['align_orth_1-96']==pos]
 if (len(set(wtlist)) > 2):
#	 print wtlist
	 cmat = yyy.T.corr('pearson')
#	 print cmat
	 n,n = cmat.shape
	 tril_values = cmat.as_matrix()[np.tril_indices(n,-1)]
	 avg_corr = tril_values.mean()
	 #print avg_corr
	 avg_corrs.append((pos,avg_corr,n,0.5*n*(n-1)))
	 #print aa,avg_corr
	 all_corrs.append(tril_values)

all_corrs = np.hstack(all_corrs)

avg_corrs_df = pd.DataFrame(avg_corrs,columns=['pos','r','num','num_triangle'])
avg_corrs_df.to_csv('./numbers_for_paper/avg_corrs_ident_pos_different_wt_%s.csv' %(version),index=False)
# # # avg_corrs_df
print "Median corr(R) over all identical structural positions, different AA: ", np.median(all_corrs)
print "Mean corr(R) over all identical structural positions, different AA: ", all_corrs.mean()
print "Mean corr(R) over all identical structural positions, different AA (not weighted): ", avg_corrs_df['r'].mean()

with open('./numbers_for_paper/avg_corrs_ident_pos_different_wt_%s.csv' %(version),'a') as fp:
 fp.write('\nmedian weighted, %.4f\nmean weighted, %.4f\nmean non-weighted, %.4f'%(np.median(all_corrs), all_corrs.mean(),avg_corrs_df['r'].mean()))


cors_dict['aligned_pos_diff_wt_aa'] = pd.Series(all_corrs)


# ##########################################
# # corrs averaging by structurally identical positions bab 22 ...
# ###########################################

avg_corrs = []
all_corrs = []
for pos in range(1,25):
 yyy = dat[aacids][dat['align_bab_1-24']==pos]
 cmat = yyy.T.corr('pearson')
 n,n = cmat.shape
 tril_values = cmat.as_matrix()[np.tril_indices(n,-1)]
 avg_corr = tril_values.mean()
 avg_corrs.append((pos,avg_corr,n,0.5*n*(n-1)))
 #print aa,avg_corr
 all_corrs.append(tril_values)

all_corrs = np.hstack(all_corrs)

avg_corrs_df = pd.DataFrame(avg_corrs,columns=['pos','r','num','num_triangle'])
avg_corrs_df.to_csv('./numbers_for_paper/avg_corrs_ident_pos22_%s.csv' %(version),index=False)
avg_corrs_df
print "Median corr(R) over all identical structural positions: ", np.median(all_corrs)
print "Mean corr(R) over all identical structural positions: ", all_corrs.mean()
print "Mean corr(R) over all identical structural positions(not weighted): ", avg_corrs_df['r'].mean()

with open('./numbers_for_paper/avg_corrs_ident_pos22_%s.csv' %(version),'a') as fp:
 fp.write('\nmedian weighted, %.4f\nmean weighted, %.4f\nmean non-weighted, %.4f'%(np.median(all_corrs), all_corrs.mean(),avg_corrs_df['r'].mean()))

cors_dict['aligned_pos_fourfold'] = pd.Series(all_corrs)




# ##########################################
# # corrs averaging by structurally identical positions and identical wt ...
# ###########################################

avg_corrs = []
all_corrs = []
for pos in range(1,97):
    grp = dat[dat['align_orth_1-96']==pos].groupby('wtaa')
    for aa,idxs in grp.groups.iteritems():
       if len(idxs)>1:
           yyy = grp.get_group(aa)[aacids]
           cmat = yyy.T.corr('pearson')
           n,n = cmat.shape
           tril_values = cmat.as_matrix()[np.tril_indices(n,-1)]
           avg_corr = tril_values.mean()
           avg_corrs.append((pos,avg_corr,n,0.5*n*(n-1)))
#          # print aa,avg_corr
           all_corrs.append(tril_values)
# #
all_corrs = np.hstack(all_corrs)

avg_corrs_df = pd.DataFrame(avg_corrs,columns=['pos','r','num','num_triangle'])
avg_corrs_df.to_csv('./numbers_for_paper/avg_corrs_ident_pos_ident_aa_%s.csv' %(version),index=False)
# # # avg_corrs_df
print "Median corr(R) over all identical structural positions and AA: ", np.median(all_corrs)
print "Mean corr(R) over all identical structural positions and AA: ", all_corrs.mean()
print "Mean corr(R) over all identical structural positions and AA(not weighted): ", avg_corrs_df['r'].mean()

with open('./numbers_for_paper/avg_corrs_ident_pos_ident_aa_%s.csv' %(version),'a') as fp:
 fp.write('\nmedian weighted, %.4f\nmean weighted, %.4f\nmean non-weighted, %.4f'%(np.median(all_corrs), all_corrs.mean(),avg_corrs_df['r'].mean()))

cors_dict['aligned_pos_identical_wt_aa'] = pd.Series(all_corrs)

# ##########################################
# # corrs averaging by structurally identical positions and similar biochemical wt ...
# ###########################################

avg_corrs = []
all_corrs = []
for pos in range(1,97):
    grp = dat[dat['align_orth_1-96']==pos].groupby('biochem')
    for feat,idxs in grp.groups.iteritems():
       if len(idxs)>1:
           yyy = grp.get_group(feat)[aacids]
           cmat = yyy.T.corr('pearson')
           n,n = cmat.shape
           tril_values = cmat.as_matrix()[np.tril_indices(n,-1)]
           avg_corr = tril_values.mean()
           avg_corrs.append((pos, feat, avg_corr,n,0.5*n*(n-1)))
#          # print feat,avg_corr
           all_corrs.append(tril_values)
# #
all_corrs = np.hstack(all_corrs)

avg_corrs_df = pd.DataFrame(avg_corrs,columns=['pos','feat','r','num','num_triangle'])
avg_corrs_df.to_csv('./numbers_for_paper/avg_corrs_ident_pos_similar_biochem_aa_%s.csv' %(version),index=False)
# # # avg_corrs_df
print "Median corr(R) over all identical structural positions and similar biochem AA: ", np.median(all_corrs)
print "Mean corr(R) over all identical structural positions and biochem similar AA: ", all_corrs.mean()
print "Mean corr(R) over all identical structural positions and biochem similar AA(not weighted): ", avg_corrs_df['r'].mean()

with open('./numbers_for_paper/avg_corrs_ident_pos_similar_biochem_aa_%s.csv' %(version),'a') as fp:
 fp.write('\nmedian weighted, %.4f\nmean weighted, %.4f\nmean non-weighted, %.4f'%(np.median(all_corrs), all_corrs.mean(),avg_corrs_df['r'].mean()))

cors_dict['aligned_pos_simil_wt_aa'] = pd.Series(all_corrs)


cors_df = pd.concat(cors_dict,names=['group',]).reset_index(level=0,drop=False).reset_index(drop=True).rename(columns={0:'corr'})
cors_df.to_csv('./numbers_for_paper/all_corrs_by_groups_%s.csv' %(version),index=False)


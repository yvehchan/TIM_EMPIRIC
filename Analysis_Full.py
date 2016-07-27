# coding: utf-8
import EMPIRIC_analysis_EXP as emp
import MSA_analysis_EXP as msa
import pandas  as pd
from scipy import stats as st
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqUtils

aacids = list(SeqUtils.IUPAC.protein.letters)


# SOME INPUT FILES WITH THE DESCRIPTION OF STRUCTURAL FEATURES,
# EMPIRIC FITNESS READOUTS AND ALIGNMENTS ...
# structural features of library positions ...
feat_fname = './Inputs/features-original.csv'
# empiric fitness readouts ...
data_fname = './Inputs/db-fitness.csv'
# (IGPS-MSA) PFam MSA ...
msa_fname = "./PFam/PF00218_full-Aln_fixed.fasta"
# (TIM-MSA)  SPalign structural alignments for tmpid = Ss,Tt,Tm ...
mssa_fname = lambda tmpid: './SPalignResults/%s_tim_aln.fasta'%tmpid
# feat_fname = '../AccessDB-files/features-original.csv'
# msa_fname  = '../AccessDB-files/uclust/PRALINE/alignment432_clust.fasta_ali'
# msa_fname = "../AccessDB-files/IGPS_seq_Praline_align/new-PRALINE-aln.fasta"
# msa_fname = "./PFam/Thermo_aln_104.fasta"

# ANALYSE IGPS-MSA ...
# do PCA and such,
mmm,msa_ws,msa_mat = msa.MSA_pipeline(feat_fname,msa_fname,dataset=['Ss','Tm','Tt'],verbose_out=True)
# msa_mat - is just a dict of 80X20 AA freq matrices from MSAs,
# but they're already log-transformed, to extract counts do:
# np.exp(msa_mat['Ss'])-1

# ANALYSE TIM-MSA ...
# do PCA and such,
sss,mssa_ws,mssa_mat = msa.MSSA_pipeline(feat_fname,msa_fname,dataset=['Ss','Tm','Tt'],verbose_out=True)


# MSA AMINO-ACID COUNTS OUTPUT ...
for tmpid in ['Ss','Tm','Tt']:
    # save raw AA counts matrices from the IGPS-MSA
    (np.exp(msa_mat[tmpid])-1).to_csv('%s_MSA_mat.csv'%tmpid)
    # save raw AA counts matrices from the TIM-MSA
    (np.exp(mssa_mat[tmpid])-1).to_csv('%s_TIM_MSA_mat.csv'%tmpid)



# ANALYSE EMPIRIC FITNESS READOUTS (PCA and such) ...
# separately at first, i.e. 3 PCAs for 80X20 matrices ...
rrrs,proj_ws = emp.EMPIRIC_pipeline_separate(feat_fname, data_fname, dataset=['Ss','Tm','Tt'])
# then jointly, i.e. 1 PCA for 240X20 matrix ...
rrrj,proj_wj = emp.EMPIRIC_pipeline_joint(feat_fname, data_fname, dataset=['Ss','Tm','Tt'])




################# DATA IS LOADED AND ANALYZED, VISUALIZATION STARTS ######################################



# EIGENVECTOR COLLINEARITY ANALYSIS ...
# ## Where all these eigenvectors are pointing anyways? ...
# 
# 'Tm' fails to be collinear with others 3rd component direction.
# For the rest of them, first 2 eigenvectors are collinear between 'Tt','Ss' and even the joint analysis.
# 
# ### EMPIRIC (joint, separate, etc.)
dot_prod_mat = np.dot(proj_ws['Ss'].T,proj_ws['Tt'])
# dot_prod_mat = np.dot(proj_wj.T,proj_ws['Ss'])
plt.imshow(dot_prod_mat,interpolation='nearest',cmap='BrBG',vmin=-0.99,vmax=0.99)
print "Upper left 4X4 piece of the dor product matrix ... "
print dot_prod_mat[:4,:4]
# np.savetxt('./extra_results_for_publication/Tt_Tm_emp_dot_prod.txt',dot_prod_mat)
plt.show()
###################################################
# ### PFam MSA 
dot_prod_mat = np.dot(msa_ws['Tm'].T,msa_ws['Tt'])
plt.imshow(dot_prod_mat,interpolation='nearest',cmap='BrBG',vmin=-0.99,vmax=0.99)
print "First 4X4 piece of the dot product matrix ..."
print dot_prod_mat[:4,:4]
np.savetxt('./extra_results_for_publication/Tt_Ss_pfam_msa_dot_prod.txt',dot_prod_mat)
###################################################
# All 3 template derived IGPS-library alignments are essentially the same ...
# 
# ### TIM MSSA ...
dot_prod_mat = np.dot(mssa_ws['Ss'].T,mssa_ws['Tt'])
plt.imshow(dot_prod_mat,interpolation='nearest',cmap='BrBG',vmin=-0.99,vmax=0.99)
print "First 4X4 piece of the dot product matrix ..."
print dot_prod_mat[:4,:4]



# # CRAZY IDEA, thermal AA usage eigenvectors are compared with the ones form EMPIRIC ...
# # ### Comparison with Thermophilic adaptation ...
# # 
# # Just to see if there are any similar sources of variance ...
# # loading projection matrices from Arch/Bact thermal adaptation analysis ...
# a1_w = np.loadtxt('../Thermal_adapt_scripts/Publication/a1w.txt')
# a2_w = np.loadtxt('../Thermal_adapt_scripts/Publication/a2w.txt')
# b1_w = np.loadtxt('../Thermal_adapt_scripts/Publication/b1w.txt')
# b2_w = np.loadtxt('../Thermal_adapt_scripts/Publication/b2w.txt')
# dot_prod_mat = np.dot(b1_w.T,b2_w)
# plt.imshow(dot_prod_mat,interpolation='nearest',cmap='BrBG',vmin=-0.99,vmax=0.99)
# print "First 4X4 piece of the dot product matrix ..."
# print dot_prod_mat[:4,:4]
# plt.show()
# # NO SENSE AT ALL, SOURCES OF VARIANCE SHARE NOTHING IN COMMON ...



#
# EXPLORATORY DATA ANALYSIS IS TO FOLLOW ...
#
######################################################
# SOME SILLY PLOTTING FUNCS, 
# could have been avoided if SEABORN was used 
######################################################
# hhh - feature columns' names and rrr is dataFrame with the data itself ...
# ## Colored scatters for EMPIRIC data
def getScatter(hhh, rrr):
    color_map_name = 'rainbow'
    #plot layout
    numfeatures = len(hhh) #skip first 22 rows: 20 aa, org-pos, and wtaa identifiers
    numColumns = 3
    numRows = ((numfeatures-1)//numColumns)+1 #ceiling function
    remainder = numfeatures - numRows*numColumns
    ##########################
    fig = plt.figure(figsize=(15, 15), facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace=.8)
    x='PC1'
    y='PC2'
    ##########################    
    for i in range(numfeatures):
        ax = fig.add_subplot(numRows,numColumns,i+1)
        cfeature = str(hhh[i])
        vmin,vmax = rrr[cfeature].min(),rrr[cfeature].max()
        ax.set_title('%s' %cfeature)
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        scat_obj = ax.scatter(rrr[x], rrr[y], s=50, cmap=color_map_name, c=rrr[cfeature],vmin=vmin,vmax=vmax)
        fig.colorbar(scat_obj)    
    # plt.suptitle('PC%d vs PC%d - %s' %(pc1_nom, pc2_nom,dataset), size = 16)
    # fig.savefig('PC%dvsPC%d-colorbyvariable-%s.pdf' %(pc1_nom, pc2_nom, dataset))
#############################################################################################################
#############################################################################################################
def sc(x,y,ax=None):
    a,b,r,p,_ = st.linregress(x,y)
    if ax is None:
        plt.plot(x,y,'ro')
        plt.plot(x,a*x+b,'b-',label="r=%.3f, p=%.3f"%(r,p))
        plt.legend(frameon=False)
    else:
        ax.plot(x,y,'ro')
        ax.plot(x,a*x+b,'b-',label="r=%.3f, p=%.3f"%(r,p))
        ax.legend(frameon=False)
#############################################################################################################
#############################################################################################################
def getLinReg(hhh, rrr, y='PC1'):
    numfeatures = len(hhh) #skip first 22 rows: 20 aa, org-pos, and wtaa identifiers
    numColumns = 3
    numRows = ((numfeatures-1)//numColumns)+1 #ceiling function
    remainder = numfeatures - numRows*numColumns
    #########################
    fig = plt.figure(figsize=(15, 15), facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace=.8)
    #########################
    for i in range(numfeatures):
        #############################################
        #
        ax = fig.add_subplot(numRows,numColumns,i+1)
        cfeature = str(hhh[i])
        ax.set_xlabel('%s' %cfeature)
        ax.set_ylabel(y)
        sc(rrr[cfeature],rrr[y],ax=ax)
        #
        ##############################################
    fig.tight_layout()


# name of the columns to use for EMPIRIC analysis
hhh_emp = ['wtKD',
           'fitness',
           'min_fitness',
           'max_fitness',
           'PC1',
           'PC2',
           'PC3',
           'PC4',
           'align_b_1-12']
# name of the columns to use for MSA analysis
hhh_msa = ['wtKD',
           'ic',
           'perc_cons',
           'gaps_profile',
           'PC1',
           'PC2',
           'PC3',
           'align_bab_1-24',
           'align_b_1-12']
# name of the columns to use for ? analysis
hhh=['align_orth_1-96',
     'align_bab_1-24',
     'align_b_1-12',
     'struct_021',
     'Loop-Type',
     'beta_inout',
     'conservation',
     'cataylsis',
     'binding',
     'ilvcluster',
     'wtKD']


# some printing ...
# make sure numberings match ...
print (rrrj['res_num'] == rrrj['pos'].apply(int)).all()


# plot some scatters ...
getScatter(hhh,rrrj); plt.show()

getScatter(hhh,mmm[mmm.organism=='Ss']); plt.show()

getScatter(hhh,sss[sss.organism=='Ss']); plt.show()


# ## EMPIRIC components are compared with features
# Compare PC from EMPIRIC PCA to see if any of the components has any clear interpretation...
getLinReg(hhh_emp, rrrs[rrrs['organism']=='Ss'], y='PC1')

getLinReg(hhh_msa, mmm[mmm['organism']=='Ss'], y='PC1')



# # SAVING RESULTS TO FILE ...
# rrrj.to_csv('EMPIRIC_All_PROCESSED.csv',index=False)
# rrrs[rrrs['organism']=='Ss'].to_csv('EMPIRIC_Ss_PROCESSED.csv',index=False)
# rrrs[rrrs['organism']=='Tm'].to_csv('EMPIRIC_Tm_PROCESSED.csv',index=False)
# rrrs[rrrs['organism']=='Tt'].to_csv('EMPIRIC_Tt_PROCESSED.csv',index=False)
# sss.to_csv('TIM_MSSA_PROCESSED.csv',index=False)
# mmm.to_csv('IGPS_PROCESSED.csv',index=False)
##################################################################################
##################################################################################
# Data saved here is used for further processing, plotting, etc...
##################################################################################
##################################################################################










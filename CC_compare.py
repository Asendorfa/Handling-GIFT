# -*- coding: utf-8 -*-

"""
Created on Wed Jun 30 10:46:36 2021

@author: Adrian Asendorf  --> adrian.asendorf@uk-koeln.de

Coded in Spyder 4. 
This is skript is specifically meant to plot already processed data. The data used
is obtained form the gifttoolbox in a mat file format. This skript is onyl suited for
this specific data structure!.
"""

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import loadmat
from collections import Counter
import os
import seaborn as sns

#==================================%MAGIC_NUMBERS%===================================================
path_dfnc_mat= 'YOUR_PROJECT_NAME_dfnc.mat'
path_post_process_mat= 'YOUR_PROJECT_NAME_dfnc_post_process.mat'
path_stats_mat='stats/YOUR_PROJECT_NAME_stats.mat'
path_subject_mat1='YOUR_PROJECT_NAME_dfnc_sub_'
path_subject_mat2='_sess_001_results.mat'
resultsdir='YOUR_OUTPUT_DIR'
grid=[0,1,3,7,9,12,21,26,36,38,44,48,50,55,57,58]
cohorts = ['GROUP1','GROUP2']
#==================================================================================================
#create a list of dictionaries
path_Cohort = [
    {   'chdir': 'PATH_TO_GIFT_PROJECT_GROUP_1',
        'group_ind': 'PATH_TO_GROUP_IND_GROUP_1'
    },
    {   'chdir': 'PATH_TO_GIFT_PROJECT_GROUP_2',
        'group_ind': 'PATH_TO_GROUP_IND_GROUP_2'
    }
]

for paths,chrt in zip(path_Cohort,cohorts):
    group_ind = pd.read_table(paths['group_ind'])
    os.chdir(paths['chdir'])
    
    #-----------------------------IMPORTING AND TRANSPOSING MAT.DATA-------------------------------------
    # 1) Data_dfnc.mat
    dfnc_mat = loadmat(path_dfnc_mat)
    dfncInfo = dfnc_mat['dfncInfo'] 
    comps    = dfncInfo['comps'][0][0][0]
    nb_comps = len(comps)
    dfnc_1   = dfncInfo[0][0][0]['comp']
    nets     = dfnc_1[0,0]
    #load the network data from matstruct
    networks=[]
    for i in range(np.shape(nets)[1]):
        networks.append(nets[0,i][0][0])
    
    
    #2) Data_postprocessing.mat
    matfile= loadmat(path_post_process_mat)
    #meta_states=matfile['meta_states_info']
    #meta_data= meta_states['tica'][0,0]
    cluster_states = matfile['clusterInfo']
    cluster_data   = cluster_states[0,0][1]
    states_info    = cluster_states[0,0][7]#
    cor_val        = cluster_data.transpose()
    fnc            = matfile['FNCcm']
    
    
    #3) Stats.mat
    matstat= loadmat(path_stats_mat)
    stat=matstat['state_vector_stats']
    frac_time_state   = stat[0,0][0]
    mean_dwell_time   = stat[0,0][1]
    transition_matrix = stat[0,0][2]
    nb_transitions    = np.concatenate(stat[0,0][3])
    
    
    #4) Import Subjects 
    all_sub_data=[]
    for sub_nb in range(1,len(states_info)+1):
        if sub_nb<10:
            path=path_subject_mat1+'00'+str(sub_nb)+path_subject_mat2
        elif sub_nb>9 and sub_nb<100:
            path=path_subject_mat1+'0'+str(sub_nb)+path_subject_mat2
        else:
            path=path_subject_mat1+str(sub_nb)+path_subject_mat2
        sub_mat=loadmat(path)
        sub_files=sub_mat['FNCdyn']
        all_sub_data.append(sub_files)
        print('->Subject '+str(sub_nb)+' sucessfully imported')
    print('-----------------------------------')
    print(f'--> Import {chrt} FINISHED')
        
    
    
    #get important numbers
    nb_clusters=len(cluster_data) 
    nb_wind=len(states_info[0][0])
    nb_tiles= len(all_sub_data[0][0])
    nb_sub= len(states_info)
    
    
    #-------------------------------------Cluster MY DATA!!!---------------------------------------------------
    #Sort the Data accoring to their windows and create correlation matrices
    
    #1)creating empty arrays for each state availabie
    for cluster in range(1,nb_clusters+1):
        vars()['state_'+str(cluster)]=[] 
    #2)sorting the windows of each subject to their respective state and calculating the median of all in windows labled state
    #representing the connectivity pattern of that state for each subject
    for sub in range(nb_sub):
        for cluster in range(1,nb_clusters+1):
            vars()['current_sub_state_'+str(cluster)]=[]
        for ind_wind,window in enumerate(states_info[sub][0]):
            for cluster in range(1,nb_clusters+1):
                if window==cluster:
                    vars()['current_sub_state_'+str(cluster)].append(all_sub_data[sub][ind_wind])
        for cluster in range(1,nb_clusters+1):
            vars()['state_'+str(cluster)].append(np.median(vars()['current_sub_state_'+str(cluster)], axis=0))
    #3)calculate the median for each state using the median correlation matrix of each PAT who was in state
    for cluster in range(1,nb_clusters+1):
        vars()[chrt+'_cleaned_state_'+str(cluster)]= [x for x in vars()['state_'+str(cluster)] if str(x) != 'nan']
        vars()[chrt+'_median_state_'+str(cluster)]=np.median(vars()[chrt+'_cleaned_state_'+str(cluster)], axis=0) 
    

def manhattan_dis (a,b):
    return sum([abs(i-j) for i,j in zip(a,b)])
def euclidian_dis(a,b):
    return np.linalg.norm(a-b)   

#write definition that transforms the 1d data in comp x comp correlation matrix
def cluster (cor_values, comps, enforce_diagonal_nan=False):
    #build a matrix of the size we want to create and fill it with ones
    nb_comps=len(comps)
    ones=np.ones((nb_comps,nb_comps))
    #make everything except the lower triangle of the matrix zero
    skeleton = np.triu(ones,1)
    #flat the matrix
    flat_skeleton = np.reshape(skeleton,(1,nb_comps*nb_comps))[0]
    #get the indexes of each one in the flattened matrix
    ones_index=[]
    for count, i in enumerate(flat_skeleton):
        if i==1:
            ones_index.append(count)
    #doing the same for transposed site
    transp_skeleton=np.tril(np.ones_like(ones,dtype=bool),-1)
    flat_transp_skeleton=np.reshape(transp_skeleton,(1,nb_comps*nb_comps))
    transp_index=[]
    for count,i in enumerate(flat_transp_skeleton[0]):
        if i==True:
            transp_index.append(count)
    for cor_val,index in zip(cor_values,ones_index):
        flat_skeleton[index]=cor_val
    mat            = np.reshape(flat_skeleton,(nb_comps,nb_comps))
    mat_transp     = np.transpose(mat)
    transp_flat    = np.reshape(mat_transp,(nb_comps*nb_comps,1))
    transp_cor_val = transp_flat[transp_flat !=0]
    #vars()[clust_name+'transp_'+str(i+1)]        = vars()[clust_name+str(i+1)]
    for corr_val2,index2 in zip(transp_cor_val,transp_index):
       flat_skeleton[index2] = corr_val2
    flat_skeleton=np.reshape(flat_skeleton,(nb_comps,nb_comps))
    np.fill_diagonal(flat_skeleton,0)
    if enforce_diagonal_nan:
        np.fill_diagonal(flat_skeleton,None)
    return flat_skeleton

#calculate a mean matrix for each network:
    
# Calculate the mean matrix for the specified rows
def calc_mean_matrix(input_matrix, grid, enforce_diagonal_zero=True, 
                     enforce_diagonal_nan=False):
    mean_matrix = []
    for i in range(len(grid[:-1]) - 1):
        start_row = grid[i]
        end_row = grid[i + 1]
        mean_row = []
        for j in range(len(grid[:-1]) - 1):
            start_col = grid[j]
            end_col = grid[j + 1]
            # Calculate mean for the specified range
            mean_val = np.nanmean(input_matrix[start_row:end_row, start_col:end_col])
            mean_row.append(mean_val)
        mean_matrix.append(mean_row)
    mean_matrix = np.array(mean_matrix)    
    if enforce_diagonal_zero:
        np.fill_diagonal(mean_matrix, 0)
    if enforce_diagonal_nan:
        np.fill_diagonal(mean_matrix, None)    
    return np.array(mean_matrix)


#-------------------------Pure tot. difference in covariance------------------------
#compare the matricies by calculating the difference in covariance
max_cor=[]
for PPMI_cluster,KFO_cluster in zip([3,4,1],[1,2,3]): #pattern of the states that are similar and 
    cluster_A = vars()['PPMI_median_state_'+str(PPMI_cluster)]
    cluster_B = vars()['KFO_median_state_'+str(KFO_cluster)]
    #calc difference
    diff_mat  = np.abs(cluster_A-cluster_B)
    #threshold
    #thresh_diff_mat = np.where(diff_mat< 0.2, 0, diff_mat)
    #cluster the data
    clustered = cluster(diff_mat,comps)
    mean_clustered = calc_mean_matrix(clustered, grid)
    #calc max cor dif
    max_cor.append(np.max(diff_mat))
    print(manhattan_dis(cluster_A, cluster_B))    
    #safe outside of the loop
    vars()['diff_mat_'+str(PPMI_cluster)+'_vs_'+str(KFO_cluster)]         = diff_mat
    vars()['diff_mat_cluster'+str(PPMI_cluster)+'_vs_'+str(KFO_cluster)]  = clustered
    vars()['mean_diff_cluster'+str(PPMI_cluster)+'_vs_'+str(KFO_cluster)] = mean_clustered

#Calculating max and min correlation vals and defining nb of clusters. 
vmax=np.max(max_cor)
    


#=================================PLOTTING=====================================
#plot the diff_ matricies
for PPMI_cluster,KFO_cluster in zip([3,4,1],[1,2,3]): #state_nbs
    clustered = vars()['diff_mat_cluster'+str(PPMI_cluster)+'_vs_'+str(KFO_cluster)]
    plt.figure(figsize=(10, 10))
    sns.heatmap(clustered, xticklabels=True, yticklabels=True, cmap="Spectral_r", 
                linewidths=0.1, linecolor='black', vmin=0, vmax=vmax, square=True, 
                cbar_kws={'label': 'Correlations (z)', 'shrink': 0.7})
    plt.title('Difference in Spearman Corr State '+str(PPMI_cluster)+' vs. State '
              +str(KFO_cluster),fontsize=25)
    [x1, x2] = plt.xlim()
    plt.hlines(grid, x1-5, x2+5, linewidths=2, color='black')
    plt.vlines(grid, *plt.ylim(), linewidths=2, color= 'black')
    
    # Adding the labels
    for pos, network in enumerate(networks):
        plt.text(-3, np.mean([grid[pos], grid[pos+1]]), network, ha='right', va='center', fontsize=12)
        plt.text(np.mean([grid[pos], grid[pos+1]]) + 0.5, nb_comps + 3, network, ha='right', va='bottom', 
                 fontsize=12, rotation=90, rotation_mode='anchor')
    #plt.show()
    plt.savefig(resultsdir+'diff_mat_cluster_'+str(PPMI_cluster)+'_vs_'+str(KFO_cluster)+'.pdf',
                bbox_inches='tight')
    plt.savefig(resultsdir+'diff_mat_cluster_'+str(PPMI_cluster)+'_vs_'+str(KFO_cluster)+'.png',
                bbox_inches='tight')
#plot mean diff mat
for PPMI_cluster,KFO_cluster in zip([3,4,1],[1,2,3]): 
    clustered = vars()['mean_diff_cluster'+str(PPMI_cluster)+'_vs_'+str(KFO_cluster)]
    plt.figure(figsize=(10, 10))
    sns.heatmap(clustered, xticklabels=True, yticklabels=True, cmap="Spectral_r", 
                linewidths=0.1, linecolor='black', vmin=0, vmax=vmax, square=True, 
                cbar_kws={'label': 'Correlations (z)', 'shrink': 0.7})
    plt.title('Mean Difference in Spearman Corr State '+str(PPMI_cluster)+' vs. State '
              +str(KFO_cluster),fontsize=25)
    [x1, x2] = plt.xlim()
    #adding labels
    for pos, network in enumerate(networks):
        plt.text(-1,pos+0.5, network, ha='right', va='center', fontsize=12)
        plt.text(pos + 0.75,len(networks)+ 0.75, network, ha='right', va='bottom', 
                 fontsize=12, rotation=90, rotation_mode='anchor')
        plt.savefig(resultsdir+'MEAN_diff_mat_cluster_'+str(PPMI_cluster)+'_vs_'+str(KFO_cluster)+'.pdf',
                    bbox_inches='tight')
        plt.savefig(resultsdir+'MEAN_diff_mat_cluster_'+str(PPMI_cluster)+'_vs_'+str(KFO_cluster)+'.png',
                    bbox_inches='tight')
    #plt.show()
  
#------------------------Statistical difference test---------------------------
#The idea of this approach is based on the Publication of Dammaraju et al. 2014
#For each tile we will compare the covariance over each tile results of each 
#individual, do a FDR correction and plot the sign vlaues.

#first we need a function to filter out the tile info of each subject:
def extract_tile_values(listofarrays, nb_tiles):
    tile_values = [[] for _ in range(nb_tiles)]
    
    # Iterate through each array
    for array in listofarrays:
        # Iterate through each tile position in the array
        for tile_index in range(nb_tiles):
            # Append the value at the current tile position to the corresponding list in tile_values
            tile_values[tile_index].append(array[tile_index])
                
    return tile_values

from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
import math

r_all = []
for PPMI_cluster,KFO_cluster in zip([3,4,1],[1,2,3]): #pattern of the states that are similar and 
    cluster_A = vars()[cohorts[0]+'_cleaned_state_'+str(PPMI_cluster)]    
    cluster_B = vars()[cohorts[1]+'_cleaned_state_'+str(KFO_cluster)]
    #get tile values of each contributing subj for the states.
    tiles_A   = extract_tile_values(cluster_A, nb_tiles)   
    tiles_B   = extract_tile_values(cluster_B, nb_tiles)
    assert len(tiles_A) == nb_tiles & len(tiles_B) == nb_tiles
    pvals = []
    rvals = []
    for tile in range(nb_tiles):
        a   = tiles_A[tile]
        b   = tiles_B[tile]
        #stat testing :MannWhiney U
        U,p = stats.mannwhitneyu(a,b)
        #calculate the Z score for the test result
        z = stats.norm.ppf(p/2) #read at Andy field Discovering Statistics using R S. 665
        #calculation on effect strength
        r= abs(z/math.sqrt(len(a)+len(b))) 
        pvals.append(p)
        rvals.append(r)
    #do an FDR correction on all those pvals
    rej, cor_pvals = fdrcorrection(pvals, alpha = 0.05, method='indep')
    #only keep rvals that are sign
    # Replace r-values corresponding to rejected null hypotheses with NaN
    for i in range(len(rvals)):
        if not rej[i]:
            rvals[i] = np.nan
    p_mat     = cluster(cor_pvals, comps)
    r_mat     = cluster(rvals, comps, enforce_diagonal_nan= True)
    #save zvals for max zval determination for scale max 
    r_all.append(rvals)
    mean_rmat = calc_mean_matrix(r_mat, grid,enforce_diagonal_nan=True)
    #safe outside of the loop
    vars()['r_mat'+str(PPMI_cluster)+'_vs_'+str(KFO_cluster)]      = r_mat
    vars()['mean_r_mat'+str(PPMI_cluster)+'_vs_'+str(KFO_cluster)] = mean_rmat
    
#Calculating max and min correlation vals and defining nb of clusters. 
vmax=np.nanmax(r_all)
vmax=1 #r is a value between 0 and 1 makes it more comaparable
    


#=================================PLOTTING=====================================
import matplotlib.colors as colors
#get only the positive half of the old colormap spectral
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

spectral = plt.get_cmap('Spectral_r')
new_spectral = truncate_colormap(spectral,0.4,1)


#plot the diff_ matricies
for PPMI_cluster,KFO_cluster in zip([3,4,1],[1,2,3]): 
    clustered = vars()['r_mat'+str(PPMI_cluster)+'_vs_'+str(KFO_cluster)]
    plt.figure(figsize=(10, 10))
    sns.heatmap(clustered, xticklabels=True, yticklabels=True, cmap=new_spectral, 
                linewidths=0.1, linecolor='black', vmin=0, vmax=vmax, square=True, 
                cbar_kws={'label': 'Effect strength (r)', 'shrink': 0.7})
    plt.title('Effect strenght Mann Whiney U Test State '+str(PPMI_cluster)+' vs. State '
              +str(KFO_cluster),fontsize=25)
    [x1, x2] = plt.xlim()
    plt.hlines(grid, x1-5, x2+5, linewidths=2, color='black')
    plt.vlines(grid, *plt.ylim(), linewidths=2, color= 'black')
    
    # Adding the labels
    for pos, network in enumerate(networks):
        plt.text(-3, np.mean([grid[pos], grid[pos+1]]), network, ha='right', va='center', fontsize=12)
        plt.text(np.mean([grid[pos], grid[pos+1]]) + 0.5, nb_comps + 3, network, ha='right', va='bottom', 
                 fontsize=12, rotation=90, rotation_mode='anchor')
    #plt.show()
    plt.savefig(resultsdir+'Effect strenght Mann Whiney U Test '+str(PPMI_cluster)+'_vs_'+str(KFO_cluster)+'.pdf',
                bbox_inches='tight')
    plt.savefig(resultsdir+'Effect strenght Mann Whiney U Test '+str(PPMI_cluster)+'_vs_'+str(KFO_cluster)+'.png',
                bbox_inches='tight')
#plot mean diff mat
for PPMI_cluster,KFO_cluster in zip([3,4,1],[1,2,3]): 
    clustered = vars()['mean_r_mat'+str(PPMI_cluster)+'_vs_'+str(KFO_cluster)]
    plt.figure(figsize=(10, 10))
    sns.heatmap(clustered, xticklabels=True, yticklabels=True, cmap=new_spectral, 
                linewidths=0.1, linecolor='black', vmin=0, vmax=vmax, square=True, 
                cbar_kws={'label': 'Mean Effect strength (r)', 'shrink': 0.7})
    plt.title('MEAN Effect strenght Mann Whiney U Test State '+str(PPMI_cluster)+' vs. State '
              +str(KFO_cluster),fontsize=25)
    [x1, x2] = plt.xlim()
    #adding labels
    for pos, network in enumerate(networks):
        plt.text(-1,pos+0.5, network, ha='right', va='center', fontsize=12)
        plt.text(pos + 0.75,len(networks)+ 0.75, network, ha='right', va='bottom', 
                 fontsize=12, rotation=90, rotation_mode='anchor')
        plt.savefig(resultsdir+'MEAN Effect strenght Mann Whiney U Test '+str(PPMI_cluster)+'_vs_'+str(KFO_cluster)+'.pdf',
                    bbox_inches='tight')
        plt.savefig(resultsdir+'MEAN Effect strenght Mann Whiney U Test'+str(PPMI_cluster)+'_vs_'+str(KFO_cluster)+'.png',
                    bbox_inches='tight')
    #plt.show()
     
        
        




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
#import seaborn as sns 



#==================================%MAGIC_NUMBERS%===================================================
os.chdir('PATH_to_Your_gift_output_dir')


path_dfnc_mat= 'YOUR_PROJECT_NAME_dfnc.mat'
path_post_process_mat= 'YOUR_PROJECT_NAME_dfnc_post_process.mat'
path_stats_mat='PATH_TO_GIFT_STATS_OUTPUT/YOUR_PROJECT_NAME_dfnc_cluster_stats.mat'
path_subject_mat1='YOUR_PROJECT_NAME_dfnc_sub_'
path_subject_mat2='_sess_001_results.mat'
resfile_path = 'SAVE_OUTPUT_OF_SKRIPT_HERE'
resultsdir=resfile_path+'plots/'
group_ind= pd.read_table('PATH_TO_GROUP_IND/group_ind.txt')#This is a text file that contains group information for each index loaded into the gifttoolbox
grid         = [0,1,3,7,9,12,21,26,36,38,44,48,50,55,57,58] #Grid koordinates for seperating Networks on Cluster
resfile_name = 'DECIDE_A_NAME_FOR_OUTPUTFILE'
demogr_path  = 'PATH_TO_YOUR_DEMOGRAPHICS_SHEET'
demogr_file  = demogr_path + 'NAME_OF_YOUR_DEMOGRAPHICS_SHEET'
used_MRI     = demogr_path + 'NAME OF SHEET CONTAINING ACTUAL IDS OF PARTICIPANTS'
#==================================================================================================



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
print('--> Import FINISHED')
    


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
    vars()['cleaned_state_'+str(cluster)]= [x for x in vars()['state_'+str(cluster)] if str(x) != 'nan']
    vars()['median_state_'+str(cluster)]=np.median(vars()['cleaned_state_'+str(cluster)], axis=0) 


#Calculating max and min correlation vals and defining nb of clusters. 
max_cor=[]
min_cor=[]
for cluster in range(1,nb_clusters+1):
    max_cor.append(np.max(vars()['median_state_'+str(cluster)]))
    min_cor.append(np.min(vars()['median_state_'+str(cluster)]))
vmax=np.max(max_cor)
vmin=np.min(min_cor)
if vmax > abs(vmin):
    v=vmax
else:
    v=vmin
clust_name='cluster_' 
print('---------------------Cluster------------------')
#print('maximal correlation values: '+str(max_cor))
print('--> Vmax= '+str(vmax) )
#print('minimal correlation values: '+str(min_cor))
print('--> Vmin= '+str(vmin) )
#print('------> V= '+str(v)

# =============================================================================
# #adjust to vmax in KFO:
# v=0.5
# 
# =============================================================================

#Calculate mean correlation vals for the absolute correlation values for 
#each state
thresh=0 #state threshold in Percent
for cluster in range(1,nb_clusters+1):
    state          = vars()['median_state_'+str(cluster)]
    filtered_state = np.array([i for i in state if i > thresh * max(state) or i < thresh * min(state)])
    absolute       = abs(filtered_state)
    mean           = round(np.mean(absolute),3)
    sumabs         = round(sum(absolute),3)
    sem            = round(stats.sem(absolute),3)
    sumnoabs       = round(sum(filtered_state),3)
    #print the results
    print('Mean absolute Correlation value state '+str(cluster)+' was '+str(mean)+
          ' (SEM='+str(sem)+') with  an absolute sum of '+str(sumabs)+' (noabs sum = '+
          str(sumnoabs)+')')

#write definition that transforms the 1d data in comp x comp correlation matrix
def cluster (cor_values, comps):
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
    return flat_skeleton


#Cluster the Corr_vectors in matrix format using the cluster function
for state in range(1,nb_clusters+1):
    vars()['cluster_'+str(state)]=pd.DataFrame(cluster(vars()['median_state_'+str(state)],comps))
    vars()['cluster_'+str(state)]=vars()['cluster_'+str(state)].set_index(comps)
    vars()['cluster_'+str(state)].columns=comps

#==========calculating the manhattan distance between cluster states:
def manhattan_dis (a,b):
    return sum([abs(i-j) for i,j in zip(a,b)])
def euclidian_dis(a,b):
    return np.linalg.norm(a-b)

for i in range(1,nb_clusters+1):
    for j in range(i + 1,nb_clusters+1):  # Start from i + 1 to avoid duplicates
    # Comparing all vectors using nested loops
        state_i=vars()['median_state_'+str(i)]
        state_j=vars()['median_state_'+str(j)]
        distance = manhattan_dis(state_i, state_j)
        print('Manhattan distance between state' + str(i) + 
              " and state" + str(j) + " is " + str(distance))

#--> calculate fnc matrix by median of all windows of each patients median 
med_wind_all_sub=[]
for sub in range(nb_sub):
    med_wind_all_sub.append(np.median(all_sub_data[sub], axis=0))

fnc_mat=pd.DataFrame(cluster(np.median(med_wind_all_sub, axis=0),comps))
fnc_mat=fnc_mat.set_index(comps)
fnc_mat.columns=comps



#-----------------------------------Grouping and mask creation----------------------------------
#creating masks for both groups

groups=groups=group_ind.columns[:-1]
for group in groups:
    vars()[str(group)+'_ind']=[int(x-1) for x in group_ind[group].dropna()]
    vars()['mask_'+str(group)]=np.zeros(nb_sub, dtype=bool)
    for i in vars()[str(group)+'_ind']:
        vars()['mask_'+str(group)][i]=True
        
    
#------------------------------------Calculate_fnc_groups-----------------------------------
for group in groups:
    median_winds=[]
    for sub in vars()[str(group)+'_ind']:
       median_winds.append(np.median(all_sub_data[sub], axis=0))
    vars()[str(group)+'_fnc_mat']=pd.DataFrame(cluster(np.median(median_winds, axis=0),comps))
    vars()[str(group)+'_fnc_mat']=vars()[str(group)+'_fnc_mat'].set_index(comps)
    vars()[str(group)+'_fnc_mat'].columns=comps
    

#-------------------------------------Dwelltime_DATA--------------------------------------------

#get info of each state and save in all states +
#count the number of states for eacht pat and calculate fractional parts
all_states_sub=[]
all_states_count=[]
states_count=[]
states_frac=[]
for sub in range(nb_sub):
     container=[]
     container2=[]
     cur_count_pd=Counter(states_info[sub][0])
     all_states_sub.append(states_info[sub][0])
     for window in range(nb_wind):
         all_states_count.append(states_info[sub][0][window])
     for states in range(1,nb_clusters+1):
         container.append(cur_count_pd[states])
         container2.append(cur_count_pd[states]/len(all_states_sub[0]))
     states_count.append(container)
     states_frac.append(container2)
all_states_count=Counter(all_states_count)
all_states_perc=[(all_states_count[x]/(nb_wind*nb_sub)) for x in range(1,nb_clusters+1)]
states_count=np.array(states_count)
states_frac=np.array(states_frac)  
states_dwell=states_count
#clacutlate the means for each group and state and store in list. 

for group in groups:
    vars()['frac_'+str(group)+'_means']=[]
    for i in range (nb_clusters):
        vars()['frac_'+str(group)+'_means'].append(np.mean(states_frac[vars()['mask_'+str(group)]][:,i]))
        
mean_frac_time=[]
for i in range (nb_clusters):
    mean_frac_time.append(np.mean(states_frac[:,i]))
    

#Count Occurances of each states for each window. 

for cluster in range(1,nb_clusters+1):
    vars()['occ_state_'+str(cluster)]=[]
    for wind in range(nb_wind):
        container=[]
        for sub in range(nb_sub):
            if states_info[sub][0][wind]==cluster:
                container.append(1)
            else:
                container.append(0)
        vars()['occ_state_'+str(cluster)].append(sum(container))

#Count Occurances of easch state for each window for each group
#Create empty lists to sort into
for cluster in range(1,nb_clusters+1):
    for group in group_ind['Group'].unique():
        vars()['occ_state_'+str(cluster)+'_'+group]=[]
#start the sorting clusterwise
for cluster in range(1,nb_clusters+1):
    for wind in range(nb_wind):
        #each window create empty container for each group
        for group in group_ind['Group'].unique():
            vars()['Container'+'_'+group]=[]
        #count the  people were in state at timepoint (window_nb)
        for sub,group in zip(range(nb_sub),group_ind['Group']):
            if states_info[sub][0][wind]==cluster:
                vars()['Container'+'_'+group].append(1)
        #add count to the empty lists 'occ_state' as percentage:
        for group in group_ind['Group'].unique():
            #define group size
            n_group=len(group_ind[group].dropna())
            #count individuals in part.s state
            n_state=sum(vars()['Container'+'_'+group])
            #calculate percentage and add to 'occ_state'
            vars()['occ_state_'+str(cluster)+'_'+group].append(
                n_state/n_group)
            


#Wrtite Function that counts the times each patient is switching beween each pair of states 
def switching_states(nb_wind,states_info,nb_clusters):
    states=range(1,nb_clusters+1)
    #define combination possibilities
    combs=[]
    for i in states:
        for j in states[i:]:
            combs.append([i,j])
    for comb in combs:
        vars()['btwn'+str(comb)+'_all']=[]
    for sub in range(len(states_info)):
        for comb in combs:
            vars()['btwn'+str(comb)]=[]
        for i in range(nb_wind-1):
            #wenn state in darauffolgenden window wechselt
            if states_info[sub][0][i]!=states_info[sub][0][i+1]:
                #zähle wenn die jewileige combination wechselt
                for comb in combs:
                    if states_info[sub][0][i]==comb[0] and states_info[sub][0][i+1]==comb[1]:
                        vars()['btwn'+str(comb)].append(1)
                    if states_info[sub][0][i]==comb[1] and states_info[sub][0][i+1]==comb[0]:
                        vars()['btwn'+str(comb)].append(1)
        for comb in combs:
            vars()['btwn'+str(comb)+'_all'].append(sum(vars()['btwn'+str(comb)]))
    switches={}
    for comb in combs:
        switches['Between_'+str(comb[0])+'n'+str(comb[1])]= vars()['btwn'+str(comb)+'_all']
   
    return switches, combs
  
#Applying function on all Data

switches_all, combs = switching_states(nb_wind,states_info,nb_clusters)
    

#switching to state
states= range(1,nb_clusters+1)
for state in states:
    vars()['switching_to_state_'+str(state)]=[]
    for sub in range(len(states_info)):
        container=[]
        for i in range(nb_wind-1):
            #For this state:For each subject if state changes from another state that is not state
            #add counter to the switching to state variable 
            if states_info[sub][0][i]!=state and states_info[sub][0][i+1]== state:
                container.append(1)
        vars()['switching_to_state_'+str(state)].append(sum(container))
        
                
            
        


#-------------------------------nb_Patients in particular state--------------------------------

for group in groups:
    vars()['sub_per_state_'+str(group)]=[]
    for i in range(nb_clusters):
        vars()['visits_state_'+str(i+1)]=[]
        vars()['sub_per_state_'+str(group)].append(len(states_dwell[vars()['mask_'+str(group)]][:,i][states_dwell[vars()['mask_'+str(group)]][:,i]!=0]))
        for sub in range(nb_sub):
            if states_dwell[:,i][sub] == 0:
                vars()['visits_state_'+str(i+1)].append(0)
            else:
                vars()['visits_state_'+str(i+1)].append(1)
for group in groups:
    vars()['sub_per_state_'+str(group)+'_perc']=[]
    for x in vars()['sub_per_state_'+str(group)]:
        vars()['sub_per_state_'+str(group)+'_perc'].append(x/len(vars()[str(group)+'_ind'])*100)



#check for monodwellers:
monodweller_equals_1=[]
for i in states_frac:
    if 1 in i:
        #print(i)
        monodweller_equals_1.append(1)
    else:
        monodweller_equals_1.append(0)

#define states monodweller stays in:
#create lists to sort 0=dweller 1=monodweller into
for i in range(nb_clusters):
    vars()['monodwells_state_'+str(i+1)]=[]
#search states_frac for 1. If 1 --> Monodweller. IF 1 not found append 0 
for i in states_frac:
    for ind,e in enumerate(i):
        if e==1:
            vars()['monodwells_state_'+str(ind+1)].append(1)
        else:
            vars()['monodwells_state_'+str(ind+1)].append(0)
    
    

monodwell_states =[]
#search states_frac for 1. If 1 --> Monodweller. IF 1 not found append 0 
for i in states_frac:
    for ind,e in enumerate(i):
        if e==1:
            monodwell_states.append(ind+1)
            break
        if ind==3 and e!=1:
           monodwell_states.append(0)


#===================================PLOTTING!================================================




import seaborn as sb
import math
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import MaxNLocator
#------------------------------Magic_plotnb------------------------------------------

#-----------------------------------------------------------------------------------
FNC_mats=[]
FRAC_times=[]

for group in groups:
    FNC_mats.append(vars()[str(group)+'_fnc_mat'])
    FRAC_times.append(vars()['frac_'+str(group)+'_means'])
    

#---------------------------------FIG1:dFC Clusters----------------------------------------


fig,axess = plt.subplots(int(math.ceil(nb_clusters/2)),2,figsize=(25,28),sharex=False)

fig.tight_layout(pad=14.0)

for i,(axes,cluster_nb) in enumerate(zip(axess.flat,range(nb_clusters))):
    sb.heatmap(vars()['cluster_'+str(i+1)], xticklabels=True, yticklabels=True,\
               ax=axes, cmap="Spectral_r",linewidths=0.1, linecolor='black', vmin=-v,vmax=v,\
                   square=True, cbar_kws={'label': 'Correlations (z)','shrink':.7})
    #cbar = axes.collections[0].colorbar
    #cbar.set_label('Correlations(z)', fontsize=13)
    #cbar.squeeze(0.7)
    axes.set_title('State '+str(i+1)+' ('+str(all_states_count[i+1])+
                   ', '+str(round(all_states_perc[i]*100,2))+'%)', fontsize=25)
    [x1,x2]=axes.get_xlim()
    axes.hlines(grid,x1-5,x2+5, linewidths=2)
    axes.vlines(grid,*axes.get_ylim(),linewidths=2)
    
    #adding the labels:
    for pos,network in enumerate(networks):
        axes.text (-3,np.mean([grid[pos],grid[pos+1]]),network, ha='right', va='center', fontsize=12)
        axes.text (np.mean([grid[pos],grid[pos+1]])+0.5,nb_comps+3,network, ha='right', va='bottom',\
                   fontsize=12, rotation=90, rotation_mode='anchor')
    
    #plt.imshow(cluster_mat_1, cmap='coolwarm', interpolation='nearest')

#Cool colors:  icefire, Spectral, YlGnBu
plt.show()
#plt.savefig('states.pdf')
fig.savefig(resultsdir+'States.pdf')
fig.savefig(resultsdir+'States.png')  


#--------------------------------FIG2:FC Cluster ------------------------------------------

fig2,axess = plt.subplots(1,1,figsize=(10,10),sharex=False)
#cbar_ax= fig.add_axes([0, 46, 51, 46])
sb.heatmap(fnc_mat, xticklabels=True, yticklabels=True,\
           cmap="Spectral_r",linewidths=0.1, linecolor='black', vmin=-v,vmax=v,\
                   square=True, cbar_kws={'label': 'Correlations (z)','shrink':.82},)
    #cbar = axes.collections[0].colorbar
    #cbar.set_label('Correlations(z)', fontsize=13)
    #cbar.squeeze(0.7)
plt.title('FNC', fontsize=25)
[x1,x2]=axes.get_xlim()
plt.hlines(grid,x1-5,x2+5, linewidths=2)
plt.vlines(grid,*axes.get_ylim(),linewidths=2)    
#adding the labels:
plt.text (-3,np.mean([grid[pos],grid[pos+1]]),network, ha='right', va='center', fontsize=12)  

for pos,network in enumerate(networks):
    plt.text (-3,np.mean([grid[pos],grid[pos+1]]),network, ha='right', va='center', fontsize=12)
    plt.text (np.mean([grid[pos],grid[pos+1]])+0.5,nb_comps+3,network, ha='right', va='bottom',\
                   fontsize=12, rotation=90, rotation_mode='anchor')
plt.show()           
fig2.savefig(resultsdir+'FNC.pdf', bbox_inches='tight') 
fig2.savefig(resultsdir+'FNC.png', bbox_inches='tight') 

#Cool colors:  icefire, Spectral, YlGnBu


#----------------------------Fig2.2 FC Cluster group-----------------------------------------------

    
for i,fnc in enumerate(FNC_mats):
    fig2,axess = plt.subplots(1,1,figsize=(10,10),sharex=False)
    #cbar_ax= fig.add_axes([0, 46, 51, 46])
    sb.heatmap(fnc, xticklabels=True, yticklabels=True,\
           cmap="Spectral_r",linewidths=0.1, linecolor='black', vmin=-v,vmax=v,\
                   square=True, cbar_kws={'label': 'Correlations (z)','shrink':.82},)
    #cbar = axes.collections[0].colorbar
    #cbar.set_label('Correlations(z)', fontsize=13)
    #cbar.squeeze(0.7)
    plt.title(str(groups[i])+' FNC', fontsize=25)
    [x1,x2]=axes.get_xlim()
    plt.hlines(grid,x1-5,x2+5, linewidths=2)
    plt.vlines(grid,*axes.get_ylim(),linewidths=2)    
    #adding the labels:
    for pos,network in enumerate(networks):
        plt.text (-3,np.mean([grid[pos],grid[pos+1]]),network, ha='right', va='center', fontsize=12)
        plt.text (np.mean([grid[pos],grid[pos+1]])+0.5,nb_comps+3,network, ha='right', va='bottom',\
                   fontsize=12, rotation=90, rotation_mode='anchor')
    plt.show()           
    fig2.savefig(resultsdir+str(groups[i])+'_FNC.pdf', bbox_inches='tight')
    fig2.savefig(resultsdir+str(groups[i])+'_FNC.png', bbox_inches='tight')



    
#---------------------------Fig5: Boxplot_Mean_Dwelltime--------------------------------

fig5,ax = plt.subplots(1,1, figsize=(8,5))

#creating list containing states for each subject
states=[]
group=[]
for i in range(nb_sub):
    states.append(list(range(1,nb_clusters+1)))
    for e in range(nb_clusters):
        group.append(group_ind['Group'][i])

#creating listto discribe group for each subject. 
dicti={'Average Dwell Time [in windows]':np.concatenate(mean_dwell_time),\
       #'Dwell Time [in windows]':np.concatenate(mean_dwell_time),\
       'states':np.concatenate(states),'Group':group}

dwelltime_df=pd.DataFrame(dicti)    
#sb.color_palette("light:#5A9", as_cmap=True)
sb.boxplot(x='states', y='Average Dwell Time [in windows]', hue='Group', data=dwelltime_df, showmeans=True,
            meanprops={"marker":"o", "markerfacecolor":"white", "markeredgecolor":"black","markersize":"5"})
plt.show()           
fig5.savefig(resultsdir+'Mean_Dwell_Time.pdf', bbox_inches='tight')
fig5.savefig(resultsdir+'Mean_Dwell_Time.png', bbox_inches='tight')


#---------------------------Fig5.5: Box_plot switches to state--------------------------------

fig5,ax = plt.subplots(1,1, figsize=(8,5))


#creating alist of all_subjects for each subject.
switching_to_state_conc=[]
for i in range(nb_sub):
    for cluster in range(nb_clusters):
        switching_to_state_conc.append(vars()['switching_to_state_'+str(cluster+1)][i])
        
    
    
dicti={'Number of Switches to state':np.array(switching_to_state_conc),\
       #'Dwell Time [in windows]':np.concatenate(mean_dwell_time),\
       'states':np.concatenate(states),'Group':np.array(group)}

dwelltime_df=pd.DataFrame(dicti)    
#sb.color_palette("light:#5A9", as_cmap=True)
sb.boxplot(x='states', y='Number of Switches to state', hue='Group', data=dwelltime_df, showmeans=True,
            meanprops={"marker":"o", "markerfacecolor":"white", "markeredgecolor":"black","markersize":"5"})
plt.show()           
fig5.savefig(resultsdir+'Switches to state.pdf', bbox_inches='tight')
fig5.savefig(resultsdir+'Switches to state.png', bbox_inches='tight')



#---------------------------Fig6: Boxplot_nb_transitions--------------------------------

fig6,ax = plt.subplots(1,1, figsize=(5,5))

#creating listto discribe group for each subject. 
dicti={'Number of Transitions':nb_transitions,'Group':group_ind['Group']}

transitions_df=pd.DataFrame(dicti)    
#sb.color_palette("light:#5A9", as_cmap=True)
#sb.violinplot(x='Group', y='Number of Transitions', data=transitions_df,\
             # split=True, palette=("YlGnBu"), cut=0)
sb.boxplot(x='Group', y='Number of Transitions', data=transitions_df, showmeans=True,\
            meanprops={"marker":"o", "markerfacecolor":"white", "markeredgecolor":"black","markersize":"6"})
plt.show()           
fig6.savefig(resultsdir+'Nb_Transitions.pdf', bbox_inches='tight') 
fig6.savefig(resultsdir+'Nb_Transitions.png', bbox_inches='tight')



#---------------------------Fig8:Barplot_ Subj_per_state_prec-----------------------------
my_colors=['mediumaquamarine','darkslategray','darkcyan']
align=['edge','center','edge']
x_fine=[0,0.25,0.5]
states=range(1,nb_clusters+1)
x_pos=[1,3,5,7,9,11,13,15,17,19,21,23]
x_pos=x_pos[0:nb_clusters]
fig8, ax= plt.subplots(1,1, figsize=(8,4))

barwidth=.5
#error=[np.std(nb_transitions[mask_hc]), np.std(nb_transitions[mask_pd])]
for i,group in enumerate(groups):
    ax.bar(x_pos,vars()['sub_per_state_'+str(group)+'_perc'], 
           align=align[i],  width=barwidth, linewidth=1.5,
           edgecolor='darkslategray', alpha=.9,color=my_colors[i])
ax.set_xlabel('state')
ax.legend(groups)
ax.set_ylabel('Subjects in state [%]')
ax.set_title('Subjects_per_state' )
ax.set_xticks(x_pos,states)
#ax.set_xlim(1,8)
plt.show()
fig8.savefig(resultsdir+'Subjects_per_state.pdf', bbox_inches='tight')
fig8.savefig(resultsdir+'Subjects_per_state.png', bbox_inches='tight')



colors=['honeydew', 'mediumaquamarine', 'lightseagreen','teal', 'steelblue','r','y','salmon','b','magenta']

#------------------Fig10:Temporal distribution of states-----------------

#for linear regression function calculating regession coefficient was written
def estimate_coef(x, y): 
	# number of observations/points 
	n = np.size(x) 

	# mean of x and y vector 
	m_x, m_y = np.mean(x), np.mean(y) 

	# calculating cross-deviation and deviation about x 
	SS_xy = np.sum(y*x) - n*m_y*m_x 
	SS_xx = np.sum(x*x) - n*m_x*m_x 

	# calculating regression coefficients 
	b_1 = SS_xy / SS_xx 
	b_0 = m_y - b_1*m_x 

	return(b_0, b_1) 

legend=[]
for i in range(1,nb_clusters+1):
    legend.append('State '+str(i))
    
fig10=plt.subplots(1,1, figsize=(10,6))
for cluster in range(1,nb_clusters+1):
    plt.scatter(range(0,nb_wind), vars()['occ_state_'+str(cluster)],\
                marker=('.'),color=colors[cluster+2], edgecolors=('grey'), linewidths=.2)
plt.legend(legend,fancybox=False)
for cluster in range(1,nb_clusters+1):    
    x=np.array(range(1,nb_wind+1))
    y=np.array(vars()['occ_state_'+str(cluster)])
    b = estimate_coef(x, y) 
    y_pred = b[0] + b[1]*x 
    plt.plot(x, y_pred, linewidth=1, color= colors[cluster+2])
    plt.xlabel('Windows')
    plt.ylabel('Number of Subjects dwelling')
    plt.text(5,np.mean(y),'y='+str(round(b[1],2))+'x+'+str(round(b[0],2)),fontsize=10\
    , alpha=1,color=colors[cluster+2])
plt.title('Temporal distribution of Number of subjects dwelling per state')      
        
plt.savefig(resultsdir+'Temporal states distribution.pdf', bbox_inches='tight')
plt.savefig(resultsdir+'Temporal states distribution.png', bbox_inches='tight')

#repeat plot for subgroups
for group in groups:    
    fig10=plt.subplots(1,1, figsize=(10,6))
    for cluster in range(1,nb_clusters+1):
        plt.scatter(range(0,nb_wind), vars()['occ_state_'+str(cluster)+'_'+group],\
                    marker=('.'),color=colors[cluster+2], edgecolors=('grey'), linewidths=.2)
    plt.legend(legend,fancybox=False)
    for cluster in range(1,nb_clusters+1):    
        x=np.array(range(1,nb_wind+1))
        y=np.array(vars()['occ_state_'+str(cluster)+'_'+group])
        b = estimate_coef(x, y) 
        y_pred = b[0] + b[1]*x 
        plt.plot(x, y_pred, linewidth=1, color= colors[cluster+2])
        plt.xlabel('Windows')
        plt.ylabel('Percentage of subjects dwelling')
        plt.text(5,np.mean(y),'y='+str(round(b[1],2))+'x+'+str(round(b[0],2)),fontsize=10\
        , alpha=1,color=colors[cluster+2])
    plt.title(group+' temporal distribution of Number of subjects dwelling per state')      
            
    plt.savefig(resultsdir+'Temporal states distribution '+group+'.pdf', bbox_inches='tight')
    plt.savefig(resultsdir+'Temporal states distribution '+group+'.png', bbox_inches='tight')



#----------------------------nb_of_switches--------------------------------------

fig12,ax = plt.subplots(1,1, figsize=(10,5))

#creating list containing states for each subject
states=[]
group=[]
for i in range(nb_sub):
    states.append(list(range(1,nb_clusters+1)))
    for e in range(nb_clusters):
        group.append(group_ind['Group'][i])

#creating listto discribe group for each subject. 
dicti={'Number of switches':[],
       #'Dwell Time [in windows]':np.concatenate(mean_dwell_time),\
       'Switching between state':[],'Group':[]}
for sub in range(nb_sub):
    for key in switches_all.keys():
        dicti['Number of switches'].append(switches_all[key][sub])
        dicti['Switching between state'].append(key[-3:])
        dicti['Group'].append(group_ind['Group'][sub])
switches_df=pd.DataFrame(dicti)    
#sb.color_palette("light:#5A9", as_cmap=True)
sb.boxplot(showmeans=True,meanline=True,meanprops={'color': 'k', 'ls': '-', 'lw': 2},medianprops={'visible': False}, 
            whiskerprops={'visible': False},zorder=10,\
            x='Switching between state', y='Number of switches', hue='Group', data=switches_df,\
            showfliers=False,showbox=False, showcaps=False,)
sb.swarmplot(x='Switching between state', y='Number of switches', hue='Group', data=switches_df,\
             linewidth=1, size=4, palette=("YlGnBu"), edgecolor=None, dodge= True )

plt.show()           
fig12.savefig(resultsdir+'nb_of_switches_btw_states.pdf', bbox_inches='tight')
fig12.savefig(resultsdir+'nb_of_switches_btw_states.png', bbox_inches='tight')



#-----------------------------------------------------------------------------------------------------
#Saving Data in one big dataFrame

dicti={}
dicti['Group']=group_ind['Group']
group_nb=[]
for i in group_ind['Group']:
    if i=='Prodromal': group_nb.append(2)
    if i=='PD': group_nb.append(3)
    if i=='Control': group_nb.append(1)
dicti['Group_nb']=group_nb
dicti['Monodwellers']=monodweller_equals_1
dicti['Monodwelled_states']= monodwell_states
for i in range(nb_clusters):
    dicti['mean_dt_st_'+str(i+1)]=mean_dwell_time[:,i]
for i in range(nb_clusters):
    dicti['monodwells_st_'+str(i+1)]=vars()['monodwells_state_'+str(i+1)]
for i in range(nb_clusters):
    dicti['switches_to_st_'+str(i+1)]=vars()['switching_to_state_'+str(i+1)]
#Add mean_dt as decile rank:
for i in range(nb_clusters):
    dicti['decile_rank_mean_dt_st_'+str(i+1)]=pd.cut(mean_dwell_time[:,i],10,labels=False, include_lowest=True)
for i in range(nb_clusters):
    dicti['frac_t_st_'+str(i+1)]=states_frac[:,i]
dicti['nb_Transitions']= nb_transitions
for i in range(1,nb_clusters+1):
    dicti['visits_state'+str(i)]=vars()['visits_state_'+str(i)]
for comb in combs:
    dicti['switch_'+str(comb[0])+'n'+str(comb[1])]= switches_all['Between_'+str(comb[0])+'n'+str(comb[1])]


#---Merge and save dfs
df_demographics = pd.read_excel(demogr_file,index_col='Subject ID')
used_IDs        = pd.read_excel(used_MRI, index_col ='Subject ID')
del used_IDs['Unnamed: 0']
#filter for IDs included in analysis
filt_df_demogra = df_demographics[df_demographics.index.isin(used_IDs.index)]
filt_df_demogra = filt_df_demogra.sort_index(axis=0, ascending= True)     
#Put all the dfsc results in one df
df_dfc          = pd.DataFrame(dicti)
#define fuction that compares lists
def lists_are_same(list1, list2):
    # Compare the lists
    return list(list1) == list(list2)
#Make shure that the sorting of both lists fits
assert lists_are_same(filt_df_demogra['Research Group'], df_dfc['Group']) == True

df_dfc          = df_dfc.set_index(filt_df_demogra.index)   
df_merged       = pd.merge(df_demographics,df_dfc,left_index=True, right_index=True)


df_merged.to_excel(resfile_path+resfile_name,index = True )  

from scipy.stats import sem
#print mean differences in days for Dat Scan for
#PD
dt_dat_PD=df_merged[df_merged['Research Group']=='PD']['Diff_days/DaT_Scan_Date'].fillna(0)
mean_dt_dat_PD= np.mean(dt_dat_PD)
std_dt_dat_PD = sem(dt_dat_PD) 
print('mean date diff in DAT dt for PD was '+str(round(mean_dt_dat_PD,1))+'_days ('+str(round(std_dt_dat_PD,2))+')')
#HC
dt_dat_HC =df_merged[df_merged['Research Group']=='Control']['Diff_days/DaT_Scan_Date'].fillna(0)
mean_dt_dat_HC=  np.mean(dt_dat_HC)
std_dt_dat_HC  = sem(dt_dat_HC) 
print('mean date diff in DAT dt for HC was '+str(round(mean_dt_dat_HC,1))+' days ('+str(round(std_dt_dat_HC,2))+')')


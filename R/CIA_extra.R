


# refine_majority_voting <- function() {
#
# }
#
#
# heatmap_group_compositions <- function() {
#
# }
#
#
#
# heatmap_grouped_distributions <- function() {
#
# }


#
#
# import numpy as np
# import pandas as pd
# import scanpy as sc
#
# def celltypist_majority_vote(data, classification_obs, groups_obs=None, min_prop=0, unassigned_label='Unassigned'):
#     """
#     A function that wraps celltypist majority vote method. It extends the the more represented cell type label (predicted by a given method) to each reference cell groups.
#     If reference cell groups are not provided it exploits scanpy.tl.leiden to over-clustering the dataset (this requires having ran scanpy.pp.neighbors or scanpy.external.pp.bbknn before).
#
#     Parameters
#     ----------
#
#
#     data: anndata.AnnData
#         an AnnData object.
#
#     classification_obs: str or list(str)
#         a string or a list of string specifying the AnnData.obs column/s where the labels assigned by the method/s of interest are stored.
#
#     groups_obs: str or None
#         a string specifying the AnnData.obs where the labels assigned by the reference method are stored. If None a over-clustering step with leiden algorithm is performed.
#
#     min_prop: float
#         for the dominant cell type within a cell group, it represent e minimum proportion of cells required to support naming of the cell group by this cell type.
#         (Default: 0, range: from 0 to 1)
#
#     unassigned_label: str
#         a string that specifies the label to assign to those cell groups in which none of the cell types reached the minimum proportion.
#
#     """
#
#     if groups_obs==None:
#         if data.n_obs < 5000:
#             resolution = 5
#         elif data.n_obs < 20000:
#             resolution = 10
#         elif data.n_obs < 40000:
#             resolution = 15
#         elif data.n_obs < 100000:
#             resolution = 20
#         elif data.n_obs < 200000:
#             resolution = 25
#         else:
#             resolution = 30
#         print('Reference annotation not selected.')
#         print('Computing over-clustering with leiden algorithm (resolution= '+str(resolution)+') ...')
#         sc.tl.leiden(data, resolution=resolution, key_added='leiden_'+str(resolution))
#         groups_obs='leiden_'+str(resolution)
#         print('Dataset has been divided into '+str(len(data.obs[groups_obs].cat.categories))+' groups accordingly with trascriptional similarities.')
#         print('')
#         print('Over-clustering result saved in AnnData.obs["'+groups_obs+'"].')
#     else:
#         print('AnnData.obs["'+groups_obs+'"] selected as reference annotation.')
#
#     print('Extending the more represented cell type label to each cell group...')
#     print('')
#     groups=np.array(data.obs[groups_obs])
#
#     if type(classification_obs)!=list:
#         classification_obs= list(classification_obs)
#     for i in classification_obs:
#         votes = pd.crosstab(data.obs[i], groups)
#         majority = votes.idxmax(axis=0)
#         freqs = (votes / votes.sum(axis=0).values).max(axis=0)
#         majority[freqs < min_prop] = 'Unassigned'
#         majority = majority[groups].reset_index()
#         majority.index =data.obs[groups_obs].index
#         majority.columns = [groups_obs, 'majority_voting']
#         majority['majority_voting'] = majority['majority_voting'].astype('category')
#         data.obs[i+' majority voting']=majority['majority_voting']
#         print('New classification labels have been stored in AnnData.obs["'+i+' majority voting"]. ')
#         print('')
#
#
# import numpy as np
# import pandas as pd
# import seaborn as sns
# import itertools
# import scipy
# import os
#
# def group_composition(data, classification_obs, groups_obs, columns_order=None, cmap='Reds', save=None):
#
#     """
#     Plots a heatmap showing the percentages of cells classified with a given method (method of interest) in cell groups defined with a different one (reference method).
#
#     Parameters
#     ----------
#
#     data: anndata.AnnData
#         an AnnData object.
#
#     classification_obs: str
#         a string specifying the AnnData.obs where the labels assigned by the method of interest are stored.
#     groups_obs: str
#         a string specifying the AnnData.obs where the labels assigned by the reference method are stored.
#     columns_order: list(str)
#         a list of strings with column names
#     cmap: str or matplotlib.colors.Colormap
#         the mapping from data values to color space.
#     save: str
#         a string specifying the file name. In the working directory, if not present, 'figures' directory will be created and a file with the prefix 'CIA_' will be saved in it.
#
#     Returns
#     -------
#
#     sns.heatmap(): AxesSubplot
#         a Axesubplot containing the percentages of cells classified with a given method in cell groups.
#     """
#     df=pd.crosstab(data.obs[groups_obs], data.obs[classification_obs])
#     df=round((df/np.array(df.sum(axis=1)).reshape(len(df.index),1))*100,2)
#     if columns_order!=None:
#         df=df[columns_order]
#         if save!=None:
#             if not os.path.exists('./figures'):
#                 os.makedirs('figures')
#             return sns.heatmap(df, cmap=cmap, annot=True).get_figure().savefig("figures/CIA_"+save)
#     return sns.heatmap(df, cmap=cmap, annot=True)
#
#
# def grouped_distributions(data, columns_obs, groups_obs, cmap='Reds', scale_medians=None, save=None):
#
#     """
#     By selecting AnnData.obs columns, this function plots a heatmap showing the medians of their values in cell groups and it prints a statistical report. For each cell group, a two-sided Wilcoxon test is perfomed to evaluate if the distribution with the highest median is different from the others. For each selected AnnData.obs columns set of values, grouped by cell groups, a two-sided Mann-Whitney U test is performed to evaluate if the distribution in the cell group having the highest median is different from the other groups distributions.
#
#     Parameters
#     ----------
#
#     data: anndata.AnnData
#         an AnnData object.
#     columns_obs: list(str)
#         a string specifying the AnnData.obs columns where the values of interest are stored.
#     groups_obs: str
#         a string specifying the AnnData.obs where the cell labels are stored.
#     cmap: str or matplotlib.colors.Colormap
#         the mapping from data values to color space.
#     scale_medians: str or None
#         a parameter to set the scaling type of median values (None, 'row-wise', 'column-wise')
#     save: str
#         a string specifying the file name. In the working directory, if not present, 'figures' directory will be created and a file with the prefix 'CIA_' will be saved in it.
#
#     Returns
#     -------
#
#     sns.heatmap(): AxesSubplot
#         a Axesubplot containing a heatmap of the score medians in cell groups.
#     """
#     grouped_df=data.obs.groupby(groups_obs).median()
#     grouped_df=grouped_df[columns_obs]
#     if scale_medians!=None:
#         if scale_medians=='row-wise':
#             grouped_df=grouped_df.transpose()/np.array(grouped_df.sum(axis=1))
#             grouped_df=grouped_df.transpose()
#         if scale_medians=='column-wise':
#             grouped_df=grouped_df/np.array(grouped_df.sum(axis=0))
#
#
#     subsets={}
#     results={}
#     print('Performing Wilcoxon test on each cell group ...')
#     combs=list(itertools.permutations(columns_obs,2))
#     count=0
#     for i in data.obs[groups_obs].cat.categories:
#         subsets[i]= data[data.obs[groups_obs]==i].obs[columns_obs]
#         pos=subsets[i].median().values.argmax()
#
#         for j in combs:
#             if ((sum(subsets[i][j[0]])!=0) & (sum(subsets[i][j[1]])!=0)):
#                 result=scipy.stats.wilcoxon(subsets[i][j[0]], subsets[i][j[1]], alternative='two-sided')
#                 if result[1] >= 0.01 and j[0]==subsets[i].median().index[pos]:
#                             count+=1
#                             print('WARNING in cell group '+i+': '+ j[0]+' values are not significantly different from '+j[1]+' values.')
#     if count==0:
#         print('For each cell group there is a distribution significantly higher than the others (p<0.01)')
#
#     print('')
#     print('Performing Mann-Whitney U test on each selected AnnData.obs column ...')
#     combs=list(itertools.permutations(data.obs[groups_obs].cat.categories,2))
#     count=0
#     for i in columns_obs:
#         sign={}
#         l=[]
#         for c in data.obs[groups_obs].cat.categories:
#             l.append(subsets[c][i].values)
#         sign[i]=pd.DataFrame(l).transpose()
#         sign[i].columns=data.obs[groups_obs].cat.categories
#         pos=sign[i].median().argmax()
#         for j in combs:
#             result=scipy.stats.mannwhitneyu(subsets[j[0]][i], subsets[j[1]][i], alternative='two-sided')
#             if result[1] >= 0.01 and j[0]==sign[i].median().index[pos]:
#                 count+=1
#                 print('WARNING in '+i+' distribution: values in '+ j[0]+' group are not significantly different from values in '+j[1]+' group')
#                 print('(p= '+str(result[1])+')')
#     if count==0:
#         print('For each distribution, there is only a cell group in which values are higher with respect to all the other groups  (p<0.01)')
#
#
#     if save!=None:
#         if not os.path.exists('./figures'):
#             os.makedirs('figures')
#         return sns.heatmap(grouped_df, cmap=cmap, annot=True).get_figure().savefig("figures/CIA_"+save)
#     return sns.heatmap(grouped_df, cmap=cmap, annot=True)

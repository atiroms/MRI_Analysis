
import numpy as np
import pandas as pd
import os
from nilearn import plotting
import nibabel as nib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.lines import Line2D


####
dir_in='424_fc_gamm_aroma_test24'
[atlas,variable,wave,model,term,sex,sign,p_threshold,id_net]=\
    ['ho112','adrenal','long','l','tanner.L',2,'both',0.005,1]

####
file_roi='D:/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/common/ROI.csv'
path_atlas='D:/NICT_WS/Dropbox/MRI_img/pnTTC/template_atlas/XCP_atlas_plot'
path_exp='D:/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP'
df_roi=pd.read_csv(file_roi)
####

df_fwep=pd.read_csv(os.path.join(path_exp,dir_in,'output','result','perm_fwep.csv'))
df_edge=pd.read_csv(os.path.join(path_exp,dir_in,'output','result','bfs_edge.csv'))

os.mkdir(os.path.join(path_exp,dir_in,'output','plot','plot_subnet'))

####
# Calculate gravity center
list_atlas=sorted(list(set(df_fwep['atlas'])))
dct_df_roi_atlas={}
for atlas in list_atlas:
    df_roi_atlas=df_roi[df_roi['atlas']==atlas].copy()
    nii_atlas=nib.load(os.path.join(path_atlas,atlas+'MNI.nii.gz'))
    [arr_gc,list_intensity]=plotting.find_parcellation_cut_coords(nii_atlas,return_label_names=True)
    for idx_row in df_roi_atlas.index:
        intensity_roi=df_roi_atlas.loc[idx_row,'intensity']
        df_roi_atlas.loc[idx_row,['gc_x','gc_y','gc_z']]=arr_gc[list_intensity.index(intensity_roi),:]
    dct_df_roi_atlas[atlas]=df_roi_atlas


# Iterate over significant networks
df_fwep_signif=df_fwep.loc[df_fwep['p_fwe']<=0.05,:]
for idx_fwep in df_fwep_signif.index:
    [atlas,variable,wave,model,term,sex,sign,p_threshold,id_net]=\
        df_fwep_signif.loc[idx_fwep,["atlas","variable","wave","model","term","sex","sign","p_threshold","id_net"]]

    # Prepare edge and node dataframe
    df_edge_subnet=df_edge.copy()
    df_edge_subnet=df_edge_subnet.loc[(df_edge['atlas']==atlas) & (df_edge['variable']==variable)\
        & (df_edge['wave']==wave) & (df_edge['model']==model) & (df_edge['term']==term)\
        & (df_edge['sex']==sex) & (df_edge['sign']==sign) & (df_edge['p_threshold']==p_threshold)\
        & (df_edge['id_net']==id_net),['from','to','estimate','F']]
    if (pd.isnull(df_edge_subnet.iloc[0]['F'])):
        df_edge_subnet=df_edge_subnet.rename(columns={'estimate':'value'})
        df_edge_subnet=df_edge_subnet.drop('F',axis=1)
        type_plot='estimate'
    else:
        df_edge_subnet=df_edge_subnet.rename(columns={'F':'value'})
        df_edge_subnet=df_edge_subnet.drop('estimate',axis=1)
        type_plot='F'
    list_node=sorted(list(set(df_edge_subnet['from'].tolist()+df_edge_subnet['to'].tolist())))
    df_roi_atlas=dct_df_roi_atlas[atlas].copy()
    df_node=df_roi_atlas.loc[df_roi_atlas['id'].isin(list_node)]
    arr_coord=df_node.loc[:,['gc_x','gc_y','gc_z']].to_numpy()

    # Prepare correlation matrix
    arr_corr=np.eye(len(list_node))
    arr_corr[arr_corr==0]=np.nan
    for idx_row in df_edge_subnet.index:
        idx_from=list_node.index(df_edge_subnet.loc[idx_row,'from'])
        idx_to=list_node.index(df_edge_subnet.loc[idx_row,'to'])
        arr_corr[idx_from,idx_to]=arr_corr[idx_to,idx_from]=df_edge_subnet.loc[idx_row,'value']

    # Plot
    arr_color=cm.nipy_spectral(np.linspace(0, 1, len(list_node)))
    if sex==1:
        label_sex='m'
    elif sex==2:
        label_sex='f'
    else:
        label_sex='mf'
    fig, axes = plt.subplots(2,1,figsize=(16,16),gridspec_kw={'height_ratios': [1,2]})                       
    display=plotting.plot_connectome(arr_corr, arr_coord,axes=axes[0],colorbar=True,
                                     edge_cmap=cm.Spectral_r,
                                     node_color=arr_color,edge_kwargs={'alpha':0.7},alpha=0.7)
    #axes[0].axis("off")
    axes[1].axis("off")
    #axes[2].axis("off")
    list_legend=[]
    for idx_node in range(len(list_node)):
        list_legend=list_legend+[Line2D([0], [0], marker='o', color='w', label=df_node.iloc[idx_node]['label'],markerfacecolor=arr_color[idx_node], markersize=5)]
    axes[1].legend(handles=list_legend, loc='upper left',bbox_to_anchor=(0, 1),ncol=4,
               fontsize='small',frameon=False)
    txt_title='atlas: '+atlas+', var: '+variable+', wave: '+str(wave)+', mod: '+model+', term: '+term+', sex: '+label_sex+', cdt: '+str(p_threshold)+', sign: '+sign+', #'+str(id_net)
    axes[0].set_title(label=txt_title,loc='center')
    fname='atl-'+atlas+'_var-'+variable+'_wav-'+str(wave)+'_mod-'+model+'_trm-'+term+'_sex-'+label_sex+'_pval-p_'+str(p_threshold)+'_sgn-'+sign+'_idx-'+str(id_net)+'_subnet.png'
    fig.savefig(os.path.join(path_exp,dir_in,'output','plot','plot_subnet',fname))
    #fig.savefig(os.path.join(path_atlas,atlas,id_roi+'.png'))

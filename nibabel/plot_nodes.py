import nibabel as nib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from nilearn.image import resample_img
from matplotlib import gridspec

####

#path_atlas='D:/atiro/Dropbox/MRI_img/pnTTC/template_atlas/XCP_atlas_plot'
#file_roi='D:/atiro/Dropbox/MRI_img/pnTTC/puberty/common/ROI.csv'
path_atlas='D:/NICT_WS/Dropbox/MRI_img/pnTTC/template_atlas/XCP_atlas_plot'
path_exp='D:/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP'
dir_in='423.2_fc_gam_cs_aroma_test3'
[atlas,variable,wave,model,term,sex,p_threshold,id_net]=\
    ['power264','corti','c1m2','l','hormone',1,0.005,1]
#['aal116','gonadal','c2m1','l','tanner',2,0.005,1]
#['power264','gonadal','c2m1','l','tanner',2,0.01,1]
#['power264','gonadal','c2m1','l','tanner',2,0.005,1]
#['power264','gonadal','c2m1','l','tanner',2,0.001,1]
#['ho112','gonadal','c2m1','l','tanner',2,0.01,1]
#['ho112','gonadal','c2m1','l','tanner',2,0.005,1]
#['ho112','gonadal','c2m1','l','tanner',2,0.001,1]

file_roi='D:/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/common/ROI.csv'
#list_atlas=['aal116','glasser360','gordon333','ho112','power264','schaefer100x7','schaefer100x17','schaefer200x7','schaefer200x17','schaefer400x7','schaefer400x17','shen268']
transparency_roi=0.3

cmap=plt.get_cmap('Wistia')

####

df_roi=pd.read_csv(file_roi)
file_template=os.path.join(path_atlas,'template.nii.gz')
nii_template=nib.load(file_template)
arr_template=nii_template.get_fdata()
arr_template=arr_template/arr_template.max()

df_roi_atlas=df_roi[df_roi['atlas']==atlas]
file_atlas=os.path.join(path_atlas,atlas+'MNI.nii.gz')
nii_atlas=nib.load(file_atlas)
if nii_atlas.shape!=nii_template.shape:
    nii_atlas=resample_img(nii_atlas,target_affine=nii_template.affine,target_shape=nii_template.shape,interpolation='nearest')
arr_atlas=nii_atlas.get_fdata()

####

df_node=pd.read_csv(os.path.join(path_exp,dir_in,'output','result','bfs_node.csv'))
df_node_subnet=df_node.copy()
df_node_subnet=df_node_subnet[(df_node['atlas']==atlas) & (df_node['variable']==variable)\
    & (df_node['wave']==wave) & (df_node['model']==model) & (df_node['term']==term)\
    & (df_node['sex']==sex) & (df_node['p_threshold']==p_threshold)\
    & (df_node['id_net']==id_net)]
df_node_subnet['degree_norm']=df_node_subnet['degree']/max(df_node_subnet['degree'])

####

arr_mask=np.zeros(shape=arr_template.shape)
arr_colored_node=np.zeros(shape=arr_template.shape+tuple([3]))

for id_row in range(df_node_subnet.shape[0]):
    [id_roi,degree_norm]=df_node_subnet.iloc[id_row][['node','degree_norm']]
    intensity_roi=int(df_roi[df_roi['id']==id_roi]['intensity'])

    arr_roi=np.zeros(shape=arr_template.shape)
    arr_roi[arr_atlas==intensity_roi]=1
    arr_roi[arr_atlas!=intensity_roi]=0

    arr_mask=arr_mask+arr_roi
    arr_colored_node[arr_atlas==intensity_roi,:]=np.asarray(cmap(degree_norm)[0:3])

    #gc_roi=np.average(np.asarray(np.where(arr_roi==1)),axis=1)
    #gc_roi_round=np.round(gc_roi)

arr_plot=np.copy(arr_template)
arr_plot[arr_mask==1]=arr_plot[arr_mask==1]*transparency_roi
arr_plot=np.stack([arr_plot,arr_plot,arr_plot],3)
arr_plot=arr_plot+arr_colored_node*(1-transparency_roi)


gs = gridspec.GridSpec(5, 8,wspace=0.05,hspace=0.05)
#gs = gridspec.GridSpec(5, 8)

#gs.update(left=0.1,right=0.1,top=0.5,bottom=0.1,wspace=0.1,hspace=0.1)
fig=plt.figure(figsize=(15,12))
idx_fig=0
for zloc in range(0,80,2):
    arr_plot_axial=np.flipud(np.fliplr(arr_plot[:,:,zloc].transpose(1,0,2)))
    idx_fig_row=idx_fig // 8
    idx_fig_col=idx_fig % 8

    plt.subplot(gs[idx_fig_row,idx_fig_col])
    plt.imshow(arr_plot_axial)
    plt.axis('off')

    idx_fig+=1

if sex==1:
    label_sex='m'
else:
    label_sex='f'
                      
plt.suptitle('atlas: '+atlas+', measure: '+variable+', wave: '+wave+', model: '+model+', expvar: '+term+', sex: '+label_sex+', p value: p<'+str(p_threshold)+', #'+str(id_net))
plt.tight_layout()

filename='atl-'+atlas+'_var-'+variable+'_wav-'+wave+'_mod-'+model+'_trm-'+term+'_sex-'+label_sex+'_pval-p_'+str(p_threshold)+'_idx-'+str(id_net)+'_subnet_node.png'
fig.savefig(os.path.join(path_exp,dir_in,'output','plot',filename))


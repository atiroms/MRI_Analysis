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
dir_in='425_fc_ca_aroma'
[atlas,variable,wave,model,term,sex,p_threshold,id_net]=\
    ['power264','corti','c1m2','l','hormone',1,0.005,1]

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

df_strength=pd.read_csv(os.path.join(path_exp,dir_in,'output','result','fc_ca_str.csv'))
df_node_subnet=df_node.copy()
df_node_subnet=df_node_subnet[(df_node['atlas']==atlas) & (df_node['variable']==variable)\
    & (df_node['wave']==wave) & (df_node['model']==model) & (df_node['term']==term)\
    & (df_node['sex']==sex) & (df_node['p_threshold']==p_threshold)\
    & (df_node['id_net']==id_net)]
df_node_subnet['degree_norm']=df_node_subnet['degree']/max(df_node_subnet['degree'])

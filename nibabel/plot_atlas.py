import nibabel as nib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from nilearn.image import resample_img
#import math
from matplotlib import gridspec

####

path_atlas='D:/atiro/Dropbox/MRI_img/pnTTC/template_atlas/XCP_atlas_plot'
file_roi='D:/atiro/Dropbox/MRI_img/pnTTC/puberty/common/ROI.csv'
list_atlas=['aal116','glasser360','gordon333','ho112','power264','schaefer100x7','schaefer100x17','schaefer200x7','schaefer200x17','schaefer400x7','schaefer400x17','shen268']
#list_atlas=['aal116']
transparency_roi=0.5

####

df_roi=pd.read_csv(file_roi)
file_template=os.path.join(path_atlas,'template.nii.gz')
nii_template=nib.load(file_template)
arr_template=nii_template.get_fdata()
arr_template=arr_template/arr_template.max()

gs = gridspec.GridSpec(2, 2, \
    width_ratios=[arr_template.shape[0],arr_template.shape[1]],\
    height_ratios=[arr_template.shape[2],arr_template.shape[1]])

for atlas in list_atlas:
    df_roi_atlas=df_roi[df_roi['atlas']==atlas]
    file_atlas=os.path.join(path_atlas,atlas+'MNI.nii.gz')
    nii_atlas=nib.load(file_atlas)
    if nii_atlas.shape!=nii_template.shape:
        nii_atlas=resample_img(nii_atlas,target_affine=nii_template.affine,target_shape=nii_template.shape,interpolation='nearest')
    arr_atlas=nii_atlas.get_fdata()
    for id_row in range(df_roi_atlas.shape[0]):
        [id_roi,label_roi,intensity_roi]=df_roi_atlas.iloc[id_row][['id','label','intensity']]
        arr_roi=np.copy(arr_atlas)
        arr_roi[arr_atlas==intensity_roi]=1
        arr_roi[arr_atlas!=intensity_roi]=0

        gc_roi=np.average(np.asarray(np.where(arr_roi==1)),axis=1)
        gc_roi_round=np.round(gc_roi)

        arr_plot=np.copy(arr_template)
        arr_plot[arr_roi==1]=arr_plot[arr_roi==1]*transparency_roi
        arr_plot=np.stack([arr_plot,arr_plot,arr_plot],3)
        arr_plot[:,:,:,0]=arr_plot[:,:,:,0]+arr_roi*(1-transparency_roi)

        arr_plot_saggital=np.flipud(arr_plot[int(gc_roi_round[0]),:,:].transpose(1,0,2))
        arr_plot_coronal=np.fliplr(np.flipud(arr_plot[:,int(gc_roi_round[1]),:].transpose(1,0,2)))
        arr_plot_axial=np.flipud(np.fliplr(arr_plot[:,:,int(gc_roi_round[2])].transpose(1,0,2)))

        fig=plt.figure(figsize=(5,5.5))
        plt.subplot(gs[0,1])
        plt.imshow(arr_plot_saggital)
        plt.axis('off')
        plt.subplot(gs[0,0])
        plt.imshow(arr_plot_coronal)
        plt.axis('off')
        plt.subplot(gs[1,0])
        plt.imshow(arr_plot_axial)
        plt.axis('off')
        plt.suptitle(id_roi+'\n'+label_roi)
        plt.tight_layout()
        fig.savefig(os.path.join(path_atlas,atlas,id_roi+'.png'))

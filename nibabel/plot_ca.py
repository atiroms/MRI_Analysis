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
path_atlas='D:/atiro/Dropbox/MRI_img/pnTTC/template_atlas/XCP_atlas_plot'
path_exp='D:/atiro/Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP'
dir_in='425_fc_ca_aroma'

file_roi='D:/atiro/Dropbox/MRI_img/pnTTC/puberty/common/ROI.csv'
#list_atlas=['aal116','glasser360','gordon333','ho112','power264','schaefer100x7','schaefer100x17','schaefer200x7','schaefer200x17','schaefer400x7','schaefer400x17','shen268']
list_atlas=['aal116','gordon333','ho112','power264','schaefer100x17','schaefer200x17','schaefer400x17','shen268']
transparency_roi=0.3

#cmap=plt.get_cmap('Wistia')
#cmap=plt.get_cmap('viridis')
#cmap=plt.get_cmap('plasma')
#cmap=plt.get_cmap('cividis')
cmap=plt.get_cmap('cool')

####

df_roi=pd.read_csv(file_roi)
file_template=os.path.join(path_atlas,'template.nii.gz')
nii_template=nib.load(file_template)
arr_template=nii_template.get_fdata()
arr_template=arr_template/arr_template.max()

#df_strength=pd.read_csv(os.path.join(path_exp,dir_in,'output','result','fc_ca_str.csv'))
df_strength=pd.read_csv(os.path.join(path_exp,dir_in,'output','result','fc_ca_str.csv'),converters={'wave_mri':str})

[atlas,wave_mri, method, sex, dim,comp]=['ho112','1','pca','female',40,2]

for atlas in list_atlas:
    df_roi_atlas=df_roi[df_roi['atlas']==atlas]
    file_atlas=os.path.join(path_atlas,atlas+'MNI.nii.gz')
    nii_atlas=nib.load(file_atlas)
    if nii_atlas.shape!=nii_template.shape:
        nii_atlas=resample_img(nii_atlas,target_affine=nii_template.affine,target_shape=nii_template.shape,interpolation='nearest')
    arr_atlas=nii_atlas.get_fdata()

    for wave_mri in ['1','2','2-1']:
        for method in ['pca','ica']:
            for sex in ['male','female','both']:
                for dim in [10,20,40]:
                    df_strength_subset=df_strength.copy()
                    df_strength_subset=df_strength_subset[(df_strength['atlas']==atlas)\
                    & (df_strength['wave_mri']==wave_mri) & (df_strength['method']==method)\
                    & (df_strength['sex']==sex) & (df_strength['dim']==dim)]
                    if len(df_strength_subset)>0:
                        for comp in range(1,dim+1):
                            df_strength_comp=df_strength_subset.copy()
                            df_strength_comp=df_strength_comp[df_strength_comp['comp']==comp]

                            df_strength_comp['strength_norm']=df_strength_comp['strength']/max(df_strength_comp['strength'])
                            
                            arr_mask=np.zeros(shape=arr_template.shape)
                            arr_colored_node=np.zeros(shape=arr_template.shape+tuple([3]))

                            for id_row in range(df_strength_comp.shape[0]):
                                [id_roi,strength_norm]=df_strength_comp.iloc[id_row][['node','strength_norm']]
                                intensity_roi=int(df_roi[df_roi['id']==id_roi]['intensity'])

                                arr_roi=np.zeros(shape=arr_template.shape)
                                arr_roi[arr_atlas==intensity_roi]=1
                                arr_roi[arr_atlas!=intensity_roi]=0

                                arr_mask=arr_mask+arr_roi
                                arr_colored_node[arr_atlas==intensity_roi,:]=np.asarray(cmap(strength_norm)[0:3])

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

                            plt.suptitle('atlas: '+atlas+', method: '+method+', wave: '+str(wave_mri)+', sex: '+sex+', dim: '+str(dim).zfill(3)+', comp: '+str(comp).zfill(3))

                            filename='atl-'+atlas+'_method-'+method+'_ses-'+str(wave_mri)+'_sex-'+sex+'_dim-'+str(dim).zfill(3)+'_comp-'+str(comp).zfill(3)+'_fc_ca_str.png'
                            fig.savefig(os.path.join(path_exp,dir_in,'output','plot',filename))
                                                    
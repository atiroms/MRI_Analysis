# Test nilearn
# original

import numpy as np
import pandas as pd
import os
from nilearn import plotting
from nilearn.image import resample_img
import nibabel as nib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.lines import Line2D


####

#path_atlas='D:/atiro/Dropbox/MRI_img/pnTTC/template_atlas/XCP_atlas_plot'
#file_roi='D:/atiro/Dropbox/MRI_img/pnTTC/puberty/common/ROI.csv'
path_atlas='D:/NICT_WS/Dropbox/MRI_img/pnTTC/template_atlas/XCP_atlas_plot'
path_exp='D:/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP'
#dir_in='423.2_fc_gam_cs_aroma_test4'
dir_in='424_fc_gamm_aroma_test24'
[atlas,variable,wave,model,term,sex,sign,p_threshold,id_net]=\
    ['ho112','adrenal','long','l','tanner.L',2,'both',0.005,1]

file_roi='D:/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/common/ROI.csv'
#list_atlas=['aal116','glasser360','gordon333','ho112','power264','schaefer100x7','schaefer100x17','schaefer200x7','schaefer200x17','schaefer400x7','schaefer400x17','shen268']
transparency_roi=0.3

#cmap=plt.get_cmap('Wistia')
cmap=plt.get_cmap('cool')

####

# Prepare atlas image data
df_roi=pd.read_csv(file_roi)
file_template=os.path.join(path_atlas,'template.nii.gz')
nii_template=nib.load(file_template)
arr_template=nii_template.get_fdata()
arr_template=arr_template/arr_template.max()

df_roi_atlas=df_roi[df_roi['atlas']==atlas].copy()
file_atlas=os.path.join(path_atlas,atlas+'MNI.nii.gz')
nii_atlas=nib.load(file_atlas)
#if nii_atlas.shape!=nii_template.shape:
#    nii_atlas=resample_img(nii_atlas,target_affine=nii_template.affine,target_shape=nii_template.shape,interpolation='nearest')
arr_atlas=nii_atlas.get_fdata()

# Calculate ROI gravity center
[arr_gc,list_intensity]=plotting.find_parcellation_cut_coords(nii_atlas,return_label_names=True)
for idx_row in df_roi_atlas.index:
    intensity_roi=df_roi_atlas.loc[idx_row,'intensity']
    df_roi_atlas.loc[idx_row,['gc_x','gc_y','gc_z']]=arr_gc[list_intensity.index(intensity_roi),:]


#for id_row in df_roi_atlas.index:
#    [id_roi,label_roi,intensity_roi]=df_roi_atlas.loc[id_row][['id','label','intensity']]
#    arr_roi=np.copy(arr_atlas)
#    arr_roi[arr_atlas==intensity_roi]=1
#    arr_roi[arr_atlas!=intensity_roi]=0
#    gc_roi=np.average(np.asarray(np.where(arr_roi==1)),axis=1)
#    df_roi_atlas.loc[id_row,['gc_x','gc_y','gc_z']]=gc_roi

# Prepare edge and node dataframe
df_edge=pd.read_csv(os.path.join(path_exp,dir_in,'output','result','bfs_edge.csv'))
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
df_node=df_roi_atlas.loc[df_roi_atlas['id'].isin(list_node)]

# Prepare correlation matrix
arr_corr=np.eye(len(list_node))
arr_corr[arr_corr==0]=np.nan
for idx_row in df_edge_subnet.index:
    idx_from=list_node.index(df_edge_subnet.loc[idx_row,'from'])
    idx_to=list_node.index(df_edge_subnet.loc[idx_row,'to'])
    arr_corr[idx_from,idx_to]=arr_corr[idx_to,idx_from]=df_edge_subnet.loc[idx_row,'value']

# Prepare node coordinate array
arr_coord=df_node.loc[:,['gc_x','gc_y','gc_z']].to_numpy()

# Plot
fig, axs = plt.subplots(2,1,figsize=(16,16),gridspec_kw={'height_ratios': [2, 3]})                       
display=plotting.plot_connectome(arr_corr, arr_coord,axes=axs[0],colorbar=True,
                                 edge_cmap=cm.Spectral_r,title='test',
                                 node_color=node_color,edge_kwargs={'alpha':0.7},alpha=0.7)
axs[1].axis("off")
list_legend=[]
for idx_node in range(len(list_node)):
    list_legend=list_legend+[Line2D([0], [0], marker='o', color='w', label=df_node.iloc[idx_node]['label'],markerfacecolor=node_color[idx_node], markersize=5)]
fig.legend(handles=list_legend, loc='upper left',bbox_to_anchor=(0.05, 0.5),ncol=4,
           fontsize='small',frameon=False)



# https://nilearn.github.io/auto_examples/03_connectivity/plot_atlas_comparison.html#sphx-glr-auto-examples-03-connectivity-plot-atlas-comparison-py

from nilearn import datasets

yeo = datasets.fetch_atlas_yeo_2011()
print('Yeo atlas nifti image (3D) with 17 parcels and liberal mask is located '
      'at: %s' % yeo['thick_17'])

####
data = datasets.fetch_development_fmri(n_subjects=10)

print('Functional nifti images (4D, e.g., one subject) are located at : %r'
      % data['func'][0])
print('Counfound csv files (of same subject) are located at : %r'
      % data['confounds'][0])

####
from nilearn.input_data import NiftiLabelsMasker
from nilearn.connectome import ConnectivityMeasure

# ConenctivityMeasure from Nilearn uses simple 'correlation' to compute
# connectivity matrices for all subjects in a list
connectome_measure = ConnectivityMeasure(kind='correlation')

# useful for plotting connectivity interactions on glass brain
from nilearn import plotting

# create masker to extract functional data within atlas parcels
masker = NiftiLabelsMasker(labels_img=yeo['thick_17'], standardize=True,
                           memory='nilearn_cache')

# extract time series from all subjects and concatenate them
time_series = []
for func, confounds in zip(data.func, data.confounds):
    time_series.append(masker.fit_transform(func, confounds=confounds))

# calculate correlation matrices across subjects and display
correlation_matrices = connectome_measure.fit_transform(time_series)

# Mean correlation matrix across 10 subjects can be grabbed like this,
# using connectome measure object
mean_correlation_matrix = connectome_measure.mean_

# grab center coordinates for atlas labels
coordinates = plotting.find_parcellation_cut_coords(labels_img=yeo['thick_17'])

# plot connectome with 80% edge strength in the connectivity
plotting.plot_connectome(mean_correlation_matrix, coordinates,
                         edge_threshold="80%",
                         title='Yeo Atlas 17 thick (func)')



# https://nilearn.github.io/auto_examples/03_connectivity/plot_multi_subject_connectome.html
import numpy as np

from nilearn import plotting

n_subjects = 4  # subjects to consider for group-sparse covariance (max: 40)


def plot_matrices(cov, prec, title, labels):
    """Plot covariance and precision matrices, for a given processing. """

    prec = prec.copy()  # avoid side effects

    # Put zeros on the diagonal, for graph clarity.
    size = prec.shape[0]
    prec[list(range(size)), list(range(size))] = 0
    span = max(abs(prec.min()), abs(prec.max()))

    # Display covariance matrix
    plotting.plot_matrix(cov, cmap=plotting.cm.bwr,
                         vmin=-1, vmax=1, title="%s / covariance" % title,
                         labels=labels)
    # Display precision matrix
    plotting.plot_matrix(prec, cmap=plotting.cm.bwr,
                         vmin=-span, vmax=span, title="%s / precision" % title,
                         labels=labels)


from nilearn import datasets
msdl_atlas_dataset = datasets.fetch_atlas_msdl()
rest_dataset = datasets.fetch_development_fmri(n_subjects=n_subjects)

# print basic information on the dataset
print('First subject functional nifti image (4D) is at: %s' %
      rest_dataset.func[0])  # 4D data


from nilearn import image
from nilearn import input_data

# A "memory" to avoid recomputation
from joblib import Memory
mem = Memory('nilearn_cache')

masker = input_data.NiftiMapsMasker(
    msdl_atlas_dataset.maps, resampling_target="maps", detrend=True,
    high_variance_confounds=True, low_pass=None, high_pass=0.01,
    t_r=2, standardize=True, memory='nilearn_cache', memory_level=1,
    verbose=2)
masker.fit()

subject_time_series = []
func_filenames = rest_dataset.func
confound_filenames = rest_dataset.confounds
for func_filename, confound_filename in zip(func_filenames,
                                            confound_filenames):
    print("Processing file %s" % func_filename)

    region_ts = masker.transform(func_filename,
                                 confounds=confound_filename)
    subject_time_series.append(region_ts)


from nilearn.connectome import GroupSparseCovarianceCV
gsc = GroupSparseCovarianceCV(verbose=2)
gsc.fit(subject_time_series)

try:
    from sklearn.covariance import GraphicalLassoCV
except ImportError:
    # for Scitkit-Learn < v0.20.0
    from sklearn.covariance import GraphLassoCV as GraphicalLassoCV

gl = GraphicalLassoCV(verbose=2)
gl.fit(np.concatenate(subject_time_series))
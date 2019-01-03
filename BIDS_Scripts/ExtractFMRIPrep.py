




#############
# LIBRARIES #
#############

import os
import shutil
import gzip
import csv
import pandas as pd


#####################################################
# Extract fMRIPrep-preprocessed data to use in CONN #
#####################################################
'''
#path_from='/media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/08_conn/05_syn_12dof_1ses'
#path_to='/media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/08_conn'
path_from='D:/atiroms/Dropbox/MRI/pnTTC/BIDS/test_5sub/08_conn/05_syn_12dof_1ses'
path_to='D:/atiroms/Dropbox/MRI/pnTTC/BIDS/test_5sub/08_conn'
#path_from='/media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/13_conn/10_syn_12dof_1ses'
#path_to='/media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/13_conn'
prefices_nii=['anat/','anat/','anat/','anat/','ses-01/func/']
suffices_nii=['_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz',
              '_space-MNI152NLin2009cAsym_label-GM_probseg.nii.gz',
              '_space-MNI152NLin2009cAsym_label-WM_probseg.nii.gz',
              '_space-MNI152NLin2009cAsym_label-CSF_probseg.nii.gz',
              '_ses-01_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz']
paths_to_subdir=['anat/t1w','anat/gm','anat/wm','anat/csf','func']
prefix_conf='ses-01/func/'
suffix_conf='_ses-01_task-rest_desc-confounds_regressors.tsv'
'''

class ExtractFMRIPrep():
    def __init__(self,path_from=path_from,path_to=path_to,
                 prefices_nii=prefices_nii,suffices_nii=suffices_nii,
                 prefix_conf=prefix_conf,suffix_conf=suffix_conf,
                 paths_to_subdir=paths_to_subdir):
        list_dir_all = os.listdir(os.path.join(path_from,'fmriprep'))
        list_dir_sub = [d for d in list_dir_all if (os.path.isdir(os.path.join(path_from,'fmriprep',d)) and d.startswith('sub-'))]
        print('List of subjects:')
        print(list_dir_sub)
        for d in ['anat','func','conf']:
            if not os.path.exists(os.path.join(path_to,d)):
                os.makedirs(os.path.join(path_to,d))
        for d in ['t1w','gm','wm','csf']:
            if not os.path.exists(os.path.join(path_to,'anat',d)):
                os.makedirs(os.path.join(path_to,'anat',d))
        print('Starting extraction...')

        for d in list_dir_sub:
            for (prefix,suffix,dir_to) in zip(prefices_nii,suffices_nii,paths_to_subdir):
                path_nii_from=os.path.join(path_from,'fmriprep',d,(prefix + d + suffix))
                path_nii_to=os.path.join(path_to,dir_to,(d + suffix))
                shutil.copy(path_nii_from,path_nii_to)
                path_nii_to_unzip,_ = os.path.splitext(path_nii_to)
                with gzip.open(path_nii_to, 'rb') as f_in:
                    with open(path_nii_to_unzip, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(path_nii_to)

            path_from_conf=os.path.join(path_from,'fmriprep',d,(prefix_conf + d + suffix_conf))
            path_to_conf=os.path.join(path_to,'conf',(d + suffix_conf))
            shutil.copy(path_from_conf,path_to_conf)

            print('Done ' + d + '.')

        print('All done.')


################################################
# Subset and summarize confounding factor data #
################################################
path_from='/media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/08_conn/conf'
path_to='/media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/08_conn/conf_subset'
columns=['']

class SubsetConf():
    def __init__(self, path_from=path_from, path_to=path_to, columns=columns):
        list_dir = os.listdir(path_from)
        for d in list_dir:
            path_file_from=os.path.join(path_from,d)
            path_file_to,_=os.path.splitext(d)
            path_file_to=os.path.join(path_to,path_file_to + '.txt')
            pd.read_table(path_file_from)



path_file_from='D:/atiroms/Dropbox/MRI/pnTTC/BIDS/test_5sub/08_conn/conf/sub-00014_ses-01_task-rest_desc-confounds_regressors.tsv'
table_in=pd.read_table(path_file_from)
table_in[0,5]

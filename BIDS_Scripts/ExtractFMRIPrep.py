
## used to extract nii.gz files from fmriprep-processed data

##############
# PARAMETERS #
##############
path_from='/media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/08_conn/05_syn_12dof_1ses'
path_to='/media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/08_conn'
prefix_t1w=''
suffix_t1w='_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz'
prefix_bold=''
suffix_bold='_ses-01_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'
prefix_conf=
suffix_conf=


#############
# LIBRARIES #
#############
import os
import shutil
import gzip


#############
# MAIN CODE #
#############

class ExtractFMRIPrep():
    def __init__(self,path_from=path_from,path_to=path_to,
                prefix_t1w=prefix_t1w,suffix_t1w=suffix_t1w,
                prefix_bold=prefix_bold,suffix_bold=suffix_bold):
        list_dir_all = os.listdir(os.path.join(path_from,'fmriprep'))
        list_dir_sub = {d for d in list_dir_all if (os.path.isdir(os.path.join(path_from,'fmriprep',d)) and d.startswith('sub-'))}
        print('List of subjects:')
        print(list_dir_sub)
        for d in {'anat','func'}:
            if not os.path.exists(os.path.join(path_to,d)):
                os.makedirs(os.path.join(path_to,d))
        for d in list_dir_sub:
            path_from_t1w=os.path.join(path_from,'fmriprep',d,'anat',(prefix_t1w + d + suffix_t1w))
            path_to_t1w=os.path.join(path_to,'anat',(prefix_t1w + d + suffix_t1w))
            shutil.copy(path_from_t1w,path_to_t1w)
            #file_to_t1w_nii,_=os.path.splitext(prefix_t1w + d + suffix_t1w)
            #path_to_t1w_nii = os.path.join(path_to_t1w, file_to_t1w_nii)
            path_to_t1w_nii,_ = os.path.splitext(path_to_t1w)
            with gzip.open(path_from_t1w, 'rb') as f_in:
                with open(path_to_t1w_nii, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(path_to_t1w)
            path_from_bold=os.path.join(path_from,'fmriprep',d,'ses-01','func',(prefix_bold + d + suffix_bold))
            path_to_bold=os.path.join(path_to,'func',(prefix_bold + d + suffix_bold))
            shutil.copy(path_from_bold,path_to_bold)
            path_to_bold_nii,_ = os.path.splitext(path_to_bold)
            with gzip.open(path_from_bold, 'rb') as f_in:
                with open(path_to_bold_nii, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(path_to_bold)
            print('Copied and unzipped T1W and BOLD file for ' + d + '.')
        print('Done.')

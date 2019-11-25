##################################################
# Libraries
##################################################

import os
import math
import numpy as np
import pandas as pd
import shutil
import glob
from tqdm.autonotebook import tqdm

#def _copyfileobj_patched(fsrc, fdst, length=16*1024*1024):
def _copyfileobj_patched(fsrc, fdst, length=1024*1024*1024):
    """Patches shutil method to hugely improve copy speed"""
    while 1:
        buf = fsrc.read(length)
        if not buf:
            break
        fdst.write(buf)
shutil.copyfileobj = _copyfileobj_patched

##################################################
# Reorganize BIDS to folders
##################################################

class Bids2dir():
    def __init__(self,
        list_path_src=['D:/MRI_img/pnTTC/data/37_c1_bids',
                       'D:/MRI_img/pnTTC/data/38_c2_bids'],
        path_dst='D:/MRI_img/pnTTC/data/400_niigz',
        list_copy=[['ses-01/anat','**/*_ses-01_T1w.nii.gz'],
                   ['ses-01/func','**/*_ses-01_task-rest_bold.nii.gz'],
                   ['ses-01/fmap','**/*_ses-01_fieldmap?.nii.gz'],
                   ['ses-02/anat','**/*_ses-02_T1w.nii.gz'],
                   ['ses-02/func','**/*_ses-02_task-rest_bold.nii.gz'],
                   ['ses-02/fmap','**/*_ses-02_fieldmap?.nii.gz']]):
        
        print('Starting Bids2dir()')

        # Create output folder
        print('Starting to create output folder.')
        list_path_mkdir=[]
        list_path_mkdir.append(path_dst)
        list_path_mkdir.append(os.path.join(path_dst,'output'))
        for content in list_copy:
            list_path_mkdir.append(os.path.join(path_dst,'output',content[0]))
        for p in list_path_mkdir:
            if not os.path.exists(p):
                os.makedirs(p)
        print('Finished creating output folder.')

        # Copy log folder
        print('Starting to copy log folder.')
        path_log_src=os.path.join(list_path_src[0],'log')
        path_log_dst=os.path.join(path_dst,'log')
        shutil.copytree(path_log_src,path_log_dst)
        print('Finished copying log folder.')

        # Copy .nii.gz files
        for content in list_copy:
            list_file_src=[]
            for path_src in list_path_src:
                list_file_src=list_file_src+glob.glob(path_src+'/'+content[1],
                                                      recursive=True)
            dir_dst=os.path.join(path_dst,'output',content[0])
            print('\nCopying: '+content[0])
            for file_src in tqdm(list_file_src):
                shutil.copy(file_src,dir_dst)
        print('Finished Bids2dir()')

##################################################
# Subset and unzip .nii.gz, create clinical data
##################################################

class SubsetNiigz():
    def __init__(self,
        path_src='D:/MRI_img/pnTTC/data/400_niigz',
        path_dst='D:/MRI_img/pnTTC/data/401_nii',
        file_clin='C:/Users/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/common/CSUB.csv',
        list_list_crt_subset=[[1,['W1_T1QC',1],['W1_rsfMRIexist',1]],
                              [2,['W2_T1QC',1],['W2_rsfMRIexist',1]]],
        list_subdir_dst=['ses-01/anat','ses-01/func','ses-02/anat','ses-02/func']):
        
        print('Starting SubsetNiigz().')

        # Create output folder
        print('Starting to create output folder.')
        list_path_mkdir=[]
        list_path_mkdir.append(path_dst)
        list_path_mkdir.append(os.path.join(path_dst,'output'))
        for subdir in list_subdir_dst:
            list_path_mkdir.append(os.path.join(path_dst,'output',subdir))
        for p in list_path_mkdir:
            if not os.path.exists(p):
                os.makedirs(p)
        print('Finished creating output folder.')

        # Copy log folder
        print('Starting to copy log folder.')
        path_log_src=os.path.join(path_src,'log')
        path_log_dst=os.path.join(path_dst,'log')
        shutil.copytree(path_log_src,path_log_dst)
        print('Finished copying log folder.')

        # Create longitudinal clinical data
        df_clin=pd.read_csv(file_clin,encoding = 'unicode_escape')
        df_clin_long=pd.DataFrame()
        for list_crt_subset in list_list_crt_subset:
            df_clin_ses=df_clin.copy()
            df_clin_ses.insert(0,'ses',list_crt_subset[0])
            for crt_subset in list_crt_subset[1:]:
                df_clin_ses=df_clin_ses[df_clin_ses[crt_subset[0]]==crt_subset[1]]
            df_clin_long=pd.concat([df_clin_long,df_clin_ses])
        df_clin_long=df_clin_long.reset_index()
        df_clin_long.to_csv(os.path.join(path_dst,'output','df_clin.csv'),index=False)


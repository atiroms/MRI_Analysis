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

##################################################
# Libraries
##################################################

import os
import math
import numpy as np
import pandas as pd
import shutil

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
# Pickup nii files
##################################################

class Pickup():
    def __init__(self,
        path_src='',
        path_dst='',
        list_ses=[1,2],
        list_file_id=['','']
        ):

        print('Starting Pickup()')

        # Create output folder
        print('Starting to create output folder.')
        list_path_mkdir=[]
        list_path_mkdir.append(path_dst)
        list_path_mkdir.append(os.path.join(path_dst,'output'))
        for ses in list_ses:
            list_path_mkdir.append(os.path.join(path_dst,'output','ses-'+str(ses).zfill(2)))
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

        # Create list of ID list
        print('Starting to create ID list.')
        list_list_id=[]
        for file_id in list_file_id:
            with open(os.path.join(path_dst,'log',file_id), 'r') as list_id:
                list_id=list_id.readlines()
                list_id=[int(x.strip('\n')) for x in list_id]
                list_id.sort()
            print('Number of subjects file: '+file_id+' : '+str(len(list_id)))
            list_list_id.append(list_id)
        print('Finished creating ID list.')

        # Copy nifti files
        print('Starting to copy nifti files.')
        for i in range(len(list_ses)):
            ses=list_ses[i]
            list_id=list_list_id[i]
            for id_subj in list_id:
                file_img='sub-'+str(id_subj).zfill(5)+'_ses-'+str(ses).zfill(2)+'_T1w.nii'
                path_file_src=os.path.join(path_src,'output','ses-'+str(ses).zfill(2),file_img)
                path_file_dst=os.path.join(path_dst,'output','ses-'+str(ses).zfill(2),file_img)
                shutil.copyfileobj(path_file_src,path_file_dst)
                print('Finished copying file: '+file_img)
        print('Finished copying nifti files.')

        print('Finished Pickup().')
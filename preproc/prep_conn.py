##################################################
# Libraries
##################################################

import os
import shutil
import pandas as pd
import csv
import nilearn.image as nl_image
import json
import numpy as np
import pydicom
import datetime
import gzip
import glob

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
# XCP to CONN
##################################################
# Preparation for using XCP-processed data to be used in CONN

class PrepCONN():
    def __init__(self,
        path_input='C:/Users/NICT_WS/MRI/pnTTC/Preproc/46_c2_nii_acompcor',
        path_output='C:/Users/NICT_WS/MRI/pnTTC/Preproc/50_c2_conn',
        file_id='id_W2_T1QC_T1QC_new_mild_rsfMRIexist_motionQC3.csv',
        session='ses-02'
        ):

        print('Starting PrepCONN().')

        # Create experiment folder
        print('Starting to create experiment folder.')
        list_paths_mkdir=[]
        list_paths_mkdir.append(path_output)
        list_paths_mkdir.append(os.path.join(path_output,'input'))
        for p in list_paths_mkdir:
            if not os.path.exists(p):
                os.makedirs(p)
        print('Finished creating experiment folder.')

        # Copy log file
        print('Starting to copy log folder.')
        path_log_in=os.path.join(path_input,'log')
        path_log_out=os.path.join(path_output,'log')
        shutil.copytree(path_log_in,path_log_out)
        print('Finished copying log folder.')

        print('Starting to load id file.')        
        with open(os.path.join(path_log_out,file_id), 'r') as list_id:
            list_id=list_id.readlines()
            list_id=[int(x.strip('\n')) for x in list_id]
            list_id.sort()
        print('Finished loading id file.')

        print('Starting to pick-up and unzip image data.')
        for id_subj in list_id:
            name_file_input='sub-'+str(id_subj).zfill(5)+'_img_sm6Std.nii.gz'
            path_file_input=os.path.join(path_input,'output',name_file_input)
            name_file_output=session+'_sub-'+str(id_subj).zfill(5)+'.nii'
            path_file_output=os.path.join(path_output,'input',name_file_output)
            with gzip.open(path_file_input, 'rb') as img_in:
                with open(path_file_output, 'wb') as img_out:
                    shutil.copyfileobj(img_in, img_out)
            print('Finished pick-up and unzipping for subject:'+ str(id_subj))
        print('Finished pick-up and unzipping image data.')

        print('Finished PrepCONN().')


##################################################
# CONN to R original
##################################################
# convert XCP- and CONN-processed data filenames to be used in R original dataset

class PostCONN():
    def __init__(self,
        path_src='D:/MRI/pnTTC/Preproc/50_c2_conn',
        path_dst='D:/MRI/pnTTC/Preproc/53_gamm',
        file_id_subj='id_W2_T1QC_T1QC_new_mild_rsfMRIexist_motionQC3.csv',
        session='ses-02'
        ):

        print('Starting PostCONN().')

        # Create experiment folder
        print('Starting to create experiment folder.')
        list_paths_mkdir=[]
        list_paths_mkdir.append(path_dst)
        list_paths_mkdir.append(os.path.join(path_dst,'input'))
        for p in list_paths_mkdir:
            if not os.path.exists(p):
                os.makedirs(p)
        print('Finished creating experiment folder.')

        # Copy log file
        print('Starting to copy log folder.')
        path_log_in=os.path.join(path_src,'log')
        path_log_out=os.path.join(path_dst,'log')
        shutil.copytree(path_log_in,path_log_out)
        print('Finished copying log folder.')

        # Load ID file
        print('Starting to load id file.')        
        with open(os.path.join(path_log_out,file_id_subj), 'r') as list_id:
            list_id=list_id.readlines()
            list_id=[int(x.strip('\n')) for x in list_id]
            list_id.sort()
        print('Finished loading id file.')
    
        # Create dataframe for ROI name conversion
        print('Starting to load ROI conversion file')
        path_roi_conn=os.path.join(path_src,'output','conn_project01','results',
                                   'firstlevel','ANALYSIS_01','_list_sources.txt')
        path_roi_dict=os.path.join(path_src,'log','ROI.csv')
        df_roiconvert=pd.read_csv(path_roi_conn,delimiter=' = ',header=None)
        df_roiconvert.columns=['label_source','label_conn']
        df_roidict=pd.read_csv(path_roi_dict)
        df_roidict=df_roidict.loc[:,['id','label_conn']]
        df_roi=pd.merge(df_roiconvert,df_roidict,on='label_conn')

        print('Starting to copy files.')
        for i in range(len(list_id)):
            for j in range(len(df_roi)):
                file_src='corr_Subject'+str(i+1).zfill(3)+'_Condition001_'+df_roi.loc[j,'label_source']+'.nii'
                path_file_src=os.path.join(path_src,'output','conn_project01','results',
                                           'firstlevel','ANALYSIS_01',file_src)
                file_dst=session+'_sub-'+str(list_id[i]).zfill(5)+'_seed-'+df_roi.loc[j,'id']+'_corr.nii'
                path_file_dst=os.path.join(path_dst,'input',file_dst)
                shutil.copyfile(path_file_src,path_file_dst)
                print('Copied ' + file_src + ' to ' + file_dst)
        print('Finished copying files.')

        print('Finished PostCONN().')


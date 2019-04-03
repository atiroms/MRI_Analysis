##################################################
# Libraries
##################################################

import os
import shutil
import pandas as pd
import csv
#import nilearn.image as nl_image
import json
import numpy as np
#import pydicom
#import datetime
#import gzip


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
# Extract motion parameter data
##################################################
# extract motion parameter data from fMRIPrep result (FSL mcflirt)
# calculate cumulative motion

class ExtractMotion():
    def __init__(self,
        #path_input='/media/veracrypt2/MRI/pnTTC/Preproc/test_5sub/35_fmriprep_latest_syn_templateout_2mm',
        #path_output='/media/veracrypt2/MRI/pnTTC/Preproc/test_5sub/50_motion',
        path_input='/media/veracrypt1/MRI/pnTTC/Preproc/26_2_fmriprep',
        path_output='/media/veracrypt1/MRI/pnTTC/Preproc/27_2_motion',
        ses='ses-02'
        ):

        print('Starting motion parameter extraction')

        # Create experiment folder
        print('Starting to create experiment folder.')
        list_paths_mkdir=[]
        list_paths_mkdir.append(path_output)
        list_paths_mkdir.append(os.path.join(path_output,'output'))
        for p in list_paths_mkdir:
            if not os.path.exists(p):
                os.makedirs(p)
        print('Finished creating experiment folder.')

        # Copy fmriprep log file
        print('Starting to copy fMRIPrep log folder.')
        path_log_in=os.path.join(path_input,'log')
        path_log_out=os.path.join(path_output,'log')
        shutil.copytree(path_log_in,path_log_out)
        print('Finished copying fMRIPrep log folder.')

        # read and process motion parmeter data
        print('Starting to extract motion data.')
        list_dir_all = os.listdir(os.path.join(path_input,'output','fmriprep'))
        list_sub=[int(d.replace('sub-','')) for d in list_dir_all if os.path.isdir(os.path.join(path_input,'output','fmriprep',d)) and d.startswith('sub-')]
        list_sub.sort()
        df_motion=pd.DataFrame(np.nan,
                               columns=['ID_pnTTC',
                                        'trans_x_max','rot_x_max',
                                        'trans_y_max','rot_y_max',
                                        'trans_z_max','rot_z_max',
                                        'trans_x_mean','rot_x_mean',
                                        'trans_y_mean','rot_y_mean',
                                        'trans_z_mean','rot_z_mean'],
                               index=range(max(list_sub)))
        df_motion.loc[:,'ID_pnTTC']=range(1,(max(list_sub)+1))
        for i in range(len(list_sub)):
            print('extracting from subject: ',str(list_sub[i]))
            name_sub='sub-'+str(list_sub[i]).zfill(5)
            path_file_confound=name_sub+'_'+ses+'_task-rest_desc-confounds_regressors.tsv'
            path_file_confound=os.path.join(path_input,'output','fmriprep',name_sub,ses,'func',path_file_confound)
            if os.path.exists(path_file_confound):
                df_confound=pd.read_csv(path_file_confound,delimiter='\t')
                for j in ['trans','rot']:
                    for k in ['x','y','z']:
                        colname=j+'_'+k
                        ts=df_confound.loc[:,colname]
                        df_motion.loc[df_motion.loc[:,'ID_pnTTC']==list_sub[i],colname+'_max']=max(abs(ts))
                        df_motion.loc[df_motion.loc[:,'ID_pnTTC']==list_sub[i],colname+'_mean']=np.mean(ts)
            else:
                print('Confound file does not exist for subject: '+str(list_sub[i]))
        path_file_output=os.path.join(path_output,'output','motion.tsv')
        df_motion.to_csv(path_file_output,sep='\t',index=False)
        print('Finished motion parameter extraction')


##################################################
# Space ID file
##################################################
# MRIQC data merging with CSUB file
# used to space skipped IDs.

class SpaceIDFile():
    def __init__(self,
        #path_file_input='/media/veracrypt1/MRI/pnTTC/Preproc/17_extractqc_ses1_t1exist/input/group_T1w.tsv',
        #path_file_output='/media/veracrypt1/MRI/pnTTC/Preproc/17_extractqc_ses1_t1exist/output/group_T1w_spaced.tsv',
        #path_file_input='/media/veracrypt1/MRI/pnTTC/Preproc/17_extractqc_ses1_t1exist/input/group_bold.tsv',
        #path_file_output='/media/veracrypt1/MRI/pnTTC/Preproc/17_extractqc_ses1_t1exist/output/group_bold_spaced.tsv',
        #path_file_input='D:/atiroms/Dropbox/MRI/pnTTC/Info/QC_New/w2_source/group_T1w.tsv',
        #path_file_output='D:/atiroms/Dropbox/MRI/pnTTC/Info/QC_New/w2_source/w2_group_T1w_spaced.csv',
        path_file_input='D:/atiroms/Dropbox/MRI/pnTTC/Info/QC_New/w2_source/group_bold.tsv',
        path_file_output='D:/atiroms/Dropbox/MRI/pnTTC/Info/QC_New/w2_source/group_bold_spaced.csv',
        colname_id_input='bids_name',
        prefix_id_input='sub-',
        #suffix_id_input='_ses-01_T1w'
        #suffix_id_input='_ses-02_T1w'
        #suffix_id_input='_ses-01_task-rest_bold'
        suffix_id_input='_ses-02_task-rest_bold'
        ):

        df=pd.read_csv(path_file_input, delimiter='\t')
        #col_id=df.loc[:,colname_id_input]
        df['ID_pnTTC']=df.apply(lambda x: int(x.loc[colname_id_input].replace(prefix_id_input,'').replace(suffix_id_input,'')),axis=1)
        list_id_input=list(df['ID_pnTTC'])
        list_id_output=list(range(1,max(list_id_input)+1,1))
        list_id_diff=list(set(list_id_output)-set(list_id_input))
        list_id_diff.sort()
        df_space=pd.DataFrame(np.nan,columns=df.columns,index=range(len(list_id_diff)))
        df_space['ID_pnTTC']=list_id_diff
        df=pd.concat([df,df_space])
        df=df.sort_values(by='ID_pnTTC')
        cols_df=df.columns.tolist()
        cols_df.remove('ID_pnTTC')
        cols_df=['ID_pnTTC'] +cols_df
        df=df[cols_df]
        #df.to_csv(path_file_output,sep='\t',index=False)
        df.to_csv(path_file_output,index=False)
        print('All done.')



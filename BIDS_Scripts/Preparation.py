#############
# LIBRARIES #
#############

import os
import shutil
import pandas as pd
#import csv
#import json


###################################
# PICUP AND COPY FREESURFER FILES #
###################################
class Fs2Fmriprep():
    def __init__(self):
        ############
        # Parameters
        path_file_id=''
        path_in=''
        path_out=''
        ############

        with open(path_file_id, 'r') as list_id:
            list_id=list_id.readlines()
            list_id=[x.strip('\n') for x in list_id]
        for i in list_id:
            path_folder_in=os.path.join(path_in,str(i).zfill(5))
            path_folder_out=os.path.join(path_out,'sub-'+str(i).zfill(5))
            shutil.copytree(path_folder_in,path_folder_out)
            print('Copied and renamed '+ path_folder_in + '.')
        print('All done.')


###################
# XCP COHORT FILE #
###################
# Create cohort file required for xcp.

class CreateCohortfile():
    def __init__(self):
        ############
        # Parameters
        #path_out='C:/Users/atiro/Dropbox/MRI/XCP_tutorial'
        #path_file_id='C:/Users/atiro/Dropbox/MRI/XCP_tutorial/id.txt'
        path_out='/media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/12_xcp'
        path_file_id='/media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/12_xcp/id.txt'
        ############

        with open(path_file_id, 'r') as list_id:
            list_id=list_id.readlines()
            list_id=[x.strip('\n') for x in list_id]
        output_anat=pd.DataFrame(columns=['id0','img'])
        #output_func=pd.DataFrame(columns=['id0','antsct','img'])
        output_func=pd.DataFrame(columns=['id0','img'])
        for index in list_id:
            output_anat=output_anat.append(pd.Series(['sub-'+str(index).zfill(5),
                                                      'fmriprep/sub-'+str(index).zfill(5)+'/anat/sub-'+str(index).zfill(5)+'_desc-preproc_T1w.nii.gz'],
                                                     index=output_anat.columns),
                                           ignore_index=True)
            output_func=output_func.append(pd.Series(['sub-'+str(index).zfill(5),
                                                      #'xcp_output/sub-'+str(index).zfill(5)+'/struc',
                                                      '10_remini_syn_12dof/fmriprep/sub-'+str(index).zfill(5)+'/ses-01/func/sub-'+str(index).zfill(5)+'_ses-01_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'],
                                                     index=output_func.columns),
                                           ignore_index=True)
        #output_anat.to_csv(os.path.join(path_out,'anat_cohort.csv'),index=False)
        output_func.to_csv(os.path.join(path_out,'func_cohort.csv'),index=False)
        print('All done.')


###############################
# MOVE ANAT FOLDER BEFORE XCP #
###############################
# move anat files in freesurfer output as workaround of xcp file reading error.

class MoveAnat():
    def __init__(self):
        ############
        # Parameters
        path_exp='/media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/12_xcp/10_remini_syn_12dof/fmriprep'
        ############

        list_dir_all = os.listdir(path_exp)
        list_sub=[d for d in list_dir_all if os.path.isdir(os.path.join(path_exp,d)) and d.startswith('sub-')]
        list_sub.sort()
        for sub in list_sub:
            path_from=os.path.join(path_exp,sub,'anat')
            path_to=os.path.join(path_exp,sub,'ses-01','anat')
            for f in os.listdir(path_from):
                shutil.move(os.path.join(path_from,f),path_to)
            print('Moved ' + sub + '/anat contents.')
        print('All done.')
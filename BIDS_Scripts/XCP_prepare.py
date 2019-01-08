#############
# LIBRARIES #
#############

import os
#import shutil
import pandas as pd
#import csv
#import json


###############
# COHORT FILE #
###############
#path_out='C:/Users/atiro/Dropbox/MRI/XCP_tutorial'
#path_file_id='C:/Users/atiro/Dropbox/MRI/XCP_tutorial/id.txt'
path_out='/media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/12_xcp'
path_file_id='/media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/12_xcp/id.txt'

class CreateCohortfile():
    def __init__(self,path_out=path_out,path_file_id=path_file_id):
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
                                                      '10_remini_syn_12dof/fmriprep/sub-'+str(index).zfill(5)+'/func/sub-'+str(index).zfill(5)+'_ses-01_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'],
                                                     index=output_func.columns),
                                           ignore_index=True)
        #output_anat.to_csv(os.path.join(path_out,'anat_cohort.csv'),index=False)
        output_func.to_csv(os.path.join(path_out,'func_cohort.csv'),index=False)
        print('All done.')


####################
# COPY ANAT FOLDER #
####################
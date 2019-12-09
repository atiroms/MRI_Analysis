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
# Clinical and brain volume data for SPM
##################################################
class ClinVol():
    def __init__(self,
        df_clin_long,
        path_dst,
        #file_clin_long='D:/MRI_img/pnTTC/c1c2_struc/spm/02_spm/output/pickup/df_clin.csv',
        file_vol='D:/MRI_img/pnTTC/c1c2_struc/spm/02_spm/output/tissue_volume/tissue_volumes.csv',
        #file_dst='D:/MRI_img/pnTTC/c1c2_struc/spm/02_spm/output/pickup/df_covar.csv',
        list_key_asis=['Sex'],
        list_key_combine=[['Age',['W1_Age_at_MRI','W2_Age_at_MRI']],
                          ['Testosterone',['W1_Testosterone','W2_Testosterone']],
                          ['DHEA',['W1_DHEA','W2_DHEA']],
                          ['Cortisol',['W1_Cortisol','W2_Cortisol']],
                          ['DHEAS',['W1_DHEAS','W2_DHEAS']]]):

        #df_clin_long=pd.read_csv(file_clin_long,encoding = 'unicode_escape')

        df_covar=df_clin_long.loc[:,['ses','ID_pnTTC']+list_key_asis]

        # Combine session-specific data into one column
        for key_combine in list_key_combine:
            for idx_row in range(len(df_covar)):
                ses_row=int(df_covar.loc[idx_row,'ses'])
                df_covar.loc[idx_row,key_combine[0]]=df_clin_long.loc[idx_row,key_combine[1][ses_row-1]]

        # load and calculate global brain calculation data
        df_vol=pd.read_csv(file_vol,encoding='unicode_escape')
        df_vol['TBV']=df_vol['Volume1']+df_vol['Volume2']
        df_vol['ICV']=df_vol['Volume1']+df_vol['Volume2']+df_vol['Volume3']
        df_vol['ses']=[int(path_file[-16:-14]) for path_file in df_vol.loc[:,'File']]
        df_vol['ID_pnTTC']=[int(path_file[-26:-21]) for path_file in df_vol.loc[:,'File']]

        # combine global calculation into covariates data
        df_covar['TBV']=pd.Series()
        df_covar['ICV']=pd.Series()
        for idx_row in range(len(df_covar)):
            df_covar.loc[idx_row,['TBV','ICV']]=df_vol.loc[(df_vol['ses']==df_covar.loc[idx_row,'ses']) & (df_vol['ID_pnTTC']==df_covar.loc[idx_row,'ID_pnTTC']),['TBV','ICV']].values.tolist()[0]
        df_covar=df_covar.sort_values(by=['ID_pnTTC','ses'])
        df_covar=df_covar.reset_index(drop=True)

        file_dst=path_dst+'/df_covar.csv'
        df_covar.to_csv(file_dst,index=False)


##################################################
# Pickup SPM-preprocessed nii files
##################################################

class Pickup():
    def __init__(self,
        path_src='D:/MRI_img/pnTTC/c1c2_struc/spm/02_spm/output/preproc',
        file_clin='C:/Users/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/common/CSUB.csv',
        file_vol='D:/MRI_img/pnTTC/c1c2_struc/spm/02_spm/output/tissue_volume/tissue_volumes.csv',
        prefix_file='smwc1',
        suffix_file='_T1w.nii',

        #path_dst='D:/MRI_img/pnTTC/c1c2_struc/spm/02_spm/output/pickup/T1QC_Hormone',
        #list_list_crt_subset=[[1,['W1_T1QC','==1'],['Sex','>0'],['W1_Testosterone','>0']],
        #                      [2,['W2_T1QC','==1'],['Sex','>0'],['W2_Testosterone','>0']]],
        #list_key_asis=['Sex'],
        #list_key_combine=[['Age',['W1_Age_at_MRI','W2_Age_at_MRI']],
        #                  ['Testosterone',['W1_Testosterone','W2_Testosterone']],
        #                  ['DHEA',['W1_DHEA','W2_DHEA']],
        #                  ['Cortisol',['W1_Cortisol','W2_Cortisol']],
        #                  ['DHEAS',['W1_DHEAS','W2_DHEAS']]]

        path_dst='D:/MRI_img/pnTTC/c1c2_struc/spm/02_spm/output/pickup/T1QC_Hormone_Female',
        list_list_crt_subset=[[1,['W1_T1QC','==1'],['Sex','==2'],['W1_Testosterone','>0']],
                              [2,['W2_T1QC','==1'],['Sex','==2'],['W2_Testosterone','>0']]],
        list_key_asis=['Sex'],
        list_key_combine=[['Age',['W1_Age_at_MRI','W2_Age_at_MRI']],
                          ['Testosterone',['W1_Testosterone','W2_Testosterone']],
                          ['DHEA',['W1_DHEA','W2_DHEA']],
                          ['Cortisol',['W1_Cortisol','W2_Cortisol']],
                          ['DHEAS',['W1_DHEAS','W2_DHEAS']]]
                          
        #path_dst='D:/MRI_img/pnTTC/c1c2_struc/spm/02_spm/output/pickup/T1QC_Hormone_Male',
        #list_list_crt_subset=[[1,['W1_T1QC','==1'],['Sex','==1'],['W1_Testosterone','>0']],
        #                      [2,['W2_T1QC','==1'],['Sex','==1'],['W2_Testosterone','>0']]],
        #list_key_asis=['Sex'],
        #list_key_combine=[['Age',['W1_Age_at_MRI','W2_Age_at_MRI']],
        #                  ['Testosterone',['W1_Testosterone','W2_Testosterone']],
        #                  ['DHEA',['W1_DHEA','W2_DHEA']],
        #                  ['Cortisol',['W1_Cortisol','W2_Cortisol']],
        #                  ['DHEAS',['W1_DHEAS','W2_DHEAS']]]

        #path_dst='D:/MRI_img/pnTTC/c1c2_struc/spm/02_spm/output/pickup/T1QC_TannerAdrenalFemale',
        #list_list_crt_subset=[[1,['W1_T1QC','==1'],['Sex','==2'],['W1_Tanner_Female_Pubic_Hair','>0']],
        #                      [2,['W2_T1QC','==1'],['Sex','==2'],['W2_Tanner_Female_Pubic_Hair','>0']]],
        #list_key_asis=['Sex'],
        #list_key_combine=[['Age',['W1_Age_at_MRI','W2_Age_at_MRI']],
        #                  ['Tanner_Adrenal',['W1_Tanner_Female_Pubic_Hair','W2_Tanner_Female_Pubic_Hair']]]
                          
        #path_dst='D:/MRI_img/pnTTC/c1c2_struc/spm/02_spm/output/pickup/T1QC_TannerAdrenalMale',
        #list_list_crt_subset=[[1,['W1_T1QC','==1'],['Sex','==1'],['W1_Tanner_Male_Pubic_Hair','>0']],
        #                      [2,['W2_T1QC','==1'],['Sex','==1'],['W2_Tanner_Male_Pubic_Hair','>0']]],
        #list_key_asis=['Sex'],
        #list_key_combine=[['Age',['W1_Age_at_MRI','W2_Age_at_MRI']],
        #                  ['Tanner_Adrenal',['W1_Tanner_Male_Pubic_Hair','W2_Tanner_Male_Pubic_Hair']]]
        
        #path_dst='D:/MRI_img/pnTTC/c1c2_struc/spm/02_spm/output/pickup/T1QC_TannerGonadalFemale',
        #list_list_crt_subset=[[1,['W1_T1QC','==1'],['Sex','==2'],['W1_Tanner_Female_Breast','>0']],
        #                      [2,['W2_T1QC','==1'],['Sex','==2'],['W2_Tanner_Female_Breast','>0']]],
        #list_key_asis=['Sex'],
        #list_key_combine=[['Age',['W1_Age_at_MRI','W2_Age_at_MRI']],
        #                  ['Tanner_Gonadal',['W1_Tanner_Female_Breast','W2_Tanner_Female_Breast']]]
                          
        #path_dst='D:/MRI_img/pnTTC/c1c2_struc/spm/02_spm/output/pickup/T1QC_TannerGonadalMale',
        #list_list_crt_subset=[[1,['W1_T1QC','==1'],['Sex','==1'],['W1_Tanner_Male_Genitals','>0']],
        #                      [2,['W2_T1QC','==1'],['Sex','==1'],['W2_Tanner_Male_Genitals','>0']]],
        #list_key_asis=['Sex'],
        #list_key_combine=[['Age',['W1_Age_at_MRI','W2_Age_at_MRI']],
        #                  ['Tanner_Gonadal',['W1_Tanner_Male_Genitals','W2_Tanner_Male_Genitals']]]
        ):

        print('Starting Pickup_preproc()')

        # Create longitudinal clinical data
        df_clin=pd.read_csv(file_clin,encoding = 'unicode_escape')
        df_clin_long=pd.DataFrame()
        for list_crt_subset in list_list_crt_subset:
            df_clin_ses=df_clin.copy()
            df_clin_ses.insert(0,'ses',list_crt_subset[0])
            for crt_subset in list_crt_subset[1:]:
                #df_clin_ses=df_clin_ses[df_clin_ses[crt_subset[0]]==crt_subset[1]]
                df_clin_ses=df_clin_ses[eval('df_clin_ses["'+crt_subset[0]+'"]'+crt_subset[1])]
            df_clin_long=pd.concat([df_clin_long,df_clin_ses])
        df_clin_long=df_clin_long.reset_index()
        df_clin_long.to_csv(os.path.join(path_dst,'df_clin_plan.csv'),index=False)

        # Pickup nii files
        df_clin_long_copied=pd.DataFrame(columns=df_clin_long.columns)
        list_absent=[]
        for idx_row in tqdm(range(len(df_clin_long))):
            file_src=prefix_file+'sub-'+str(df_clin_long.loc[idx_row,'ID_pnTTC']).zfill(5)+'_ses-'+str(df_clin_long.loc[idx_row,'ses']).zfill(2)+suffix_file
            path_file_src=glob.glob(path_src+'/'+file_src,recursive=True)
            if len(path_file_src)>0:
                path_file_src=path_file_src[0]
                shutil.copy(path_file_src,path_dst)
                df_clin_long_copied=df_clin_long_copied.append(df_clin_long.loc[idx_row,:])
            else:
                list_absent.append(prefix_file+'sub-'+str(df_clin_long.loc[idx_row,'ID_pnTTC']).zfill(5)+'_ses-'+str(df_clin_long.loc[idx_row,'ses']).zfill(2))+suffix_file

        if len(list_absent)>0:
            print('Absent data:')
            print(list_absent)
        df_clin_long_copied.to_csv(os.path.join(path_dst,'df_clin.csv'),index=False)

        # Create dataframe of clinical and global volume measures
        df_covar=df_clin_long_copied.loc[:,['ses','ID_pnTTC']+list_key_asis]

        # Combine session-specific data into one column
        for key_combine in list_key_combine:
            for idx_row in range(len(df_covar)):
                ses_row=int(df_covar.loc[idx_row,'ses'])
                df_covar.loc[idx_row,key_combine[0]]=df_clin_long.loc[idx_row,key_combine[1][ses_row-1]]

        # load and calculate global brain measurement data
        df_vol=pd.read_csv(file_vol,encoding='unicode_escape')
        df_vol['TBV']=df_vol['Volume1']+df_vol['Volume2']
        df_vol['ICV']=df_vol['Volume1']+df_vol['Volume2']+df_vol['Volume3']
        df_vol['ses']=[int(path_file[-16:-14]) for path_file in df_vol.loc[:,'File']]
        df_vol['ID_pnTTC']=[int(path_file[-26:-21]) for path_file in df_vol.loc[:,'File']]

        # combine global calculation into covariates data
        df_covar['TBV']=pd.Series()
        df_covar['ICV']=pd.Series()
        for idx_row in range(len(df_covar)):
            df_covar.loc[idx_row,['TBV','ICV']]=df_vol.loc[(df_vol['ses']==df_covar.loc[idx_row,'ses']) & (df_vol['ID_pnTTC']==df_covar.loc[idx_row,'ID_pnTTC']),['TBV','ICV']].values.tolist()[0]
        df_covar=df_covar.sort_values(by=['ID_pnTTC','ses'])
        df_covar=df_covar.reset_index(drop=True)

        df_covar.to_csv(path_dst+'/df_covar.csv',index=False)

        print('Finished Pickup_preproc()')


##################################################
# Pickup nii files
##################################################

class Pickup_old():
    def __init__(self,
        path_src='D:/MRI/pnTTC/c1c2_struc/spm/00_acpc',
        path_dst='D:/MRI/pnTTC/c1c2_struc/spm/01_t1qc',
        list_ses=[1,2],
        list_file_id=['id_ses-01_t1qc.csv','id_ses-02_t1qc.csv']
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
                shutil.copy(path_file_src,path_file_dst)
                print('Finished copying file: '+file_img)
        print('Finished copying nifti files.')

        print('Finished Pickup().')
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
# XCP cohort file
##################################################
# XCP prepataion
# Create cohort file required for xcp.

class CreateCohortfile():
    def __init__(self,
        n_proc=1,
        #path_out='C:/Users/atiro/Dropbox/MRI/XCP_tutorial',
        #path_file_id='C:/Users/atiro/Dropbox/MRI/XCP_tutorial/id.txt',
        #path_file_out='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/33_xcp_36p_templatein/input/func_cohort.csv',
        path_dir_out='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/33_xcp_36p_templatein/input',
        path_file_id='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/33_xcp_36p_templatein/log/id_5sub.txt',
        suffix_file='_ses-01_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz',
        #suffix_file='_ses-01_task-rest_space-T1w_desc-preproc_bold.nii.gz',
        #dir_input='10_remini_syn_12dof'
        #dir_input='16_fmriprep_newfs'
        ses='ses-01'
        ):

        print('Starting to create XCP cohort file.')
        with open(path_file_id, 'r') as list_id:
            list_id=list_id.readlines()
            list_id=[int(x.strip('\n')) for x in list_id]
            list_id.sort()

        # Make list of ID lists
        n_subj=len(list_id)
        n_subj_per_proc=int(np.ceil(n_subj/n_proc))
        n_proc_floor=n_proc*n_subj_per_proc-n_subj
        n_proc_ceil=n_proc-n_proc_floor
        print('  '+str(n_subj)+' total subs, '+str(n_proc)+' total procs, '+str(n_proc_ceil)+' procs with '
              +str(n_subj_per_proc)+' subs, '+str(n_proc_floor)+' procs with '+ str(n_subj_per_proc-1)+' subs.')
        list_list_id=[]
        for id_proc in range(n_proc):
            if id_proc<n_proc_ceil:
                list_list_id.append(list_id[(id_proc*n_subj_per_proc):((id_proc+1)*n_subj_per_proc)])
            else:
                list_list_id.append(list_id[(n_subj_per_proc*n_proc_ceil+(id_proc-n_proc_ceil)*(n_subj_per_proc-1)):
                                            (n_subj_per_proc*n_proc_ceil+(id_proc-n_proc_ceil+1)*(n_subj_per_proc-1))])
        
        #output_anat=pd.DataFrame(columns=['id0','img'])
        #output_func=pd.DataFrame(columns=['id0','antsct','img'])
        for id_proc in range(n_proc):
            output_func=pd.DataFrame(columns=['id0','img'])
            for id_subj in list_list_id[id_proc]:
                #output_anat=output_anat.append(pd.Series(['sub-'+str(index).zfill(5),
                #                                          'input/fmriprep/sub-'+str(index).zfill(5)+'/anat/sub-'+str(index).zfill(5)+'_desc-preproc_T1w.nii.gz'],
                #                                         index=output_anat.columns),
                #                               ignore_index=True)
                output_func=output_func.append(pd.Series(['sub-'+str(id_subj).zfill(5),
                                                          #'xcp_output/sub-'+str(id_subj).zfill(5)+'/struc',
                                                          'input/fmriprep/sub-'+str(id_subj).zfill(5)+'/'+ses+'/func/sub-'+str(id_subj).zfill(5)+suffix_file],
                                                         index=output_func.columns),
                                               ignore_index=True)
            #output_anat.to_csv(os.path.join(path_out,'anat_cohort.csv'),index=False)
            path_file_out=os.path.join(path_dir_out,'func_cohort_'+str(id_proc).zfill(2)+'.csv')
            output_func.to_csv(path_file_out,index=False)
        print('Finished creating XCP cohort file.')


##################################################
# Move anat folder
##################################################
# XCP preparation
# move anat files in freesurfer output as workaround of xcp file reading error.

class MoveAnat():
    def __init__(self,
        #path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/22_xcp_aroma_aromain/input/fmriprep'
        path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/33_xcp_36p_templatein/input/fmriprep',
        ses='ses-01'
        ):

        print('Starting to move fMRIPrep anat folder contents.')
        list_dir_all = os.listdir(path_exp)
        list_sub=[d for d in list_dir_all if os.path.isdir(os.path.join(path_exp,d)) and d.startswith('sub-')]
        list_sub.sort()
        for sub in list_sub:
            path_from=os.path.join(path_exp,sub,'anat')
            path_to=os.path.join(path_exp,sub,ses,'anat')
            for f in os.listdir(path_from):
                shutil.move(os.path.join(path_from,f),path_to)
            print('Moved ' + sub + '/anat contents.')
        
        print('Finished moving fMRIPrep anat folder contents.')


##################################################
# Generate multiple scripts for parallel XCP
##################################################
class XCPScript():
    def __init__(self,
        n_proc=1,
        path_exp='',
        file_design='',
        path_img_xcp='',
        script=''
        ):

        output=''
        for id_proc in range(n_proc):
            str_id_proc=str(id_proc).zfill(2)
            script_proc_head='# Process No. '+ str_id_proc + '\n'
            script_proc=script
            script_proc=script_proc.replace('{path_exp}',path_exp).replace('{file_design}',file_design).replace('{path_img_xcp}',path_img_xcp)
            script_proc=script_proc.replace('{id_proc}',str_id_proc)
            script_proc_tail='\n\n\n'
            script_proc=script_proc_head+script_proc+script_proc_tail
            output=output+script_proc

        path_file_output=os.path.join(path_exp,'log','XCP_script.sh')
        file=open(path_file_output,'w')
        file.write(output)
        file.close()       


##################################################
# XCP preparation
##################################################
# joining the above three classes and some more

class XCPPrep():
    def __init__(self,
        skip_fmriprep_copy=True,
        skip_fmriprep_moveanat=True,
        n_proc=20,
        path_fmriprep='/media/veracrypt1/MRI/pnTTC/Preproc/26_2_fmriprep',
        path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/29_2_xcp_aroma',
        #file_id='w2_id_mild_1.csv',
        file_id='w2_id_mild_2_omit328.csv',
        ses='ses-02',
        #suffix_img='_ses-01_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz',
        suffix_img='_ses-02_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz',
        #suffix_img='_ses-01_task-rest_space-T1w_desc-preproc_bold.nii.gz',
        #path_folder_design='/home/atiroms/Documents/GitHub/MRI_Analysis/Preprocessing/XCP_design/accessed_on_20190131/modified',
        path_folder_design='/home/atiroms/GitHub/MRI_Analysis/preproc/XCP_design/accessed_on_20190131/modified',
        #file_design='fc-36p_spkreg_fconly_noqcfc.dsn',
        file_design='fc-aroma_fconly_noqcfc.dsn',
        #file_design='fc-acompcor_fconly_noqcfc.dsn',
        #file_design='fc-36p_spkreg_fconly.dsn',
        #file_design='fc-aroma_fconly.dsn',
        #file_design='fc-acompcor_fconly.dsn',
        #file_design='fc-acompcor.dsn',
        #file_design='fc-acompcor_fc_roiquant.dsn',
        path_img_xcp='/data/applications/xcpEngine-070-20190130.simg',
        script='singularity run --cleanenv -B {path_exp}:${HOME}/data {path_img_xcp} -d ${HOME}/data/input/{file_design} -c ${HOME}/data/input/func_cohort_{id_proc}.csv -o ${HOME}/data/output/{id_proc} -t 1 -r ${HOME}/data'
        ):

        print('Starting XCP preparation.')

        # Create experiment folder
        print('Starting to create experiment folder.')
        list_paths_mkdir=[]
        list_paths_mkdir.append(path_exp)
        for d in ['input','output']:
            list_paths_mkdir.append(os.path.join(path_exp,d))
        for i in range(n_proc):
            list_paths_mkdir.append(os.path.join(path_exp,'output',str(i).zfill(2)))
        for p in list_paths_mkdir:
            if not os.path.exists(p):
                os.makedirs(p)
        print('Finished creating experiment folder.')

        # Copy fmriprep log file
        print('Starting to copy fMRIPrep log folder.')
        path_log_in=os.path.join(path_fmriprep,'log')
        path_log_out=os.path.join(path_exp,'log')
        shutil.copytree(path_log_in,path_log_out)
        print('Finished copying fMRIPrep log folder.')

        # Create XCP cohort file
        _=CreateCohortfile(n_proc=n_proc,
                           #path_file_out=os.path.join(path_exp,'input/func_cohort.csv'),
                           path_dir_out=os.path.join(path_exp,'input'),
                           path_file_id=os.path.join(path_exp,'log',file_id),
                           suffix_file=suffix_img,
                           ses=ses)

        # Copy XCP design file from Git local repsitory
        print('Starting to copy XCP design file.')
        path_dsn_in=os.path.join(path_folder_design,file_design)
        path_dsn_out=os.path.join(path_exp,'input',file_design)
        shutil.copy(path_dsn_in,path_dsn_out)
        print('Finished copying XCP design file.')

        # Generate XCP script
        print('Starting to generate XCP script.')
        _=XCPScript(n_proc, path_exp,file_design,path_img_xcp,script)
        print('Finished generating XCP script.')

        # Copy fMRIPrep preprocessed data to /input folder
        if not skip_fmriprep_copy:
            print('Starting to copy fMRIPrep folder.')
            path_fmriprep_in=os.path.join(path_fmriprep,'output/fmriprep')
            path_fmriprep_out=os.path.join(path_exp,'input/fmriprep')
            shutil.copytree(path_fmriprep_in,path_fmriprep_out)
            print('Finished copying fmriprep folder.')

        # Move fmriprep /anat folder contents (workaround of XCP bug)
        if not skip_fmriprep_moveanat:
            print("Starting to move contents of /anat folder.")
            _=MoveAnat(path_exp=os.path.join(path_exp,'input/fmriprep'),
                       ses=ses)
            print("Finished moving contents of /anat folder.")

        print('Finished XCP preparation.')


##################################################
# Pickup and unzip nii.gz data
##################################################
# CONN preparation

class PickupUnzip():
    def __init__(self,
        path_input='D:/atiroms/MRI/pnTTC/pnTTC1_rsfMRI_C/31_xcpout_aroma',
        path_output='D:/atiroms/MRI/pnTTC/pnTTC1_rsfMRI_C/33_conn_aroma',
        path_file_clinical='D:/atiroms/Dropbox/MRI/pnTTC/Puberty/Stats/CommonData/CSUB_W1_T1QC_new_mild_rsfMRIexist_motionQC3.csv'
        ):

        print('Starting pick-up and unzipping of .nii.gz data.')

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

        print('Starting to load clinical data.')        
        df_clinical=pd.read_csv(path_file_clinical,encoding='cp932')
        print('Finished loading clinical data.')

        print('Starting to pick-up and unzip image data.')
        for id_subj in df_clinical.loc[:,'ID_pnTTC']:
            name_file_input='sub-'+str(id_subj).zfill(5)+'_img_sm6Std.nii.gz'
            path_file_input=os.path.join(path_input,'output','norm',name_file_input)
            name_file_output='sub-'+str(id_subj).zfill(5)+'_img_sm6Std.nii'
            path_file_output=os.path.join(path_output,'input',name_file_output)
            with gzip.open(path_file_input, 'rb') as img_in:
                with open(path_file_output, 'wb') as img_out:
                    shutil.copyfileobj(img_in, img_out)
            print('Finished pick-up and unzipping for subject:'+ str(id_subj))
        print('Finished pick-up and unzipping image data.')

        print('Finished.')


##################################################
# Extract Normalized NIfTI data
##################################################
# Extract XCP-preprocessed NIfTI data

class ExtractNifti():
    def __init__(self,
        #path_input='Q:/MRI/pnTTC/Preproc/test_5sub/44_xcp_parallel',
        #path_output='Q:/MRI/pnTTC/Preproc/test_5sub/51_nifti'
        #path_input='Q:/MRI/pnTTC/Preproc/23_1_xcp_acompcor',
        #path_output='D:/atiroms/MRI/pnTTC/pnTTC1_rsfMRI_C/30_1_conn_acompcor'
        #path_input='P:/MRI/pnTTC/Preproc/23_2_xcp_acompcor',
        #path_output='D:/atiroms/MRI/pnTTC/pnTTC1_rsfMRI_C/30_2_conn_acompcor'
        #path_input='Q:/MRI/pnTTC/Preproc/22_1_xcp_aroma',
        #path_output='D:/atiroms/MRI/pnTTC/pnTTC1_rsfMRI_C/31_1_conn_aroma'
        #path_input='P:/MRI/pnTTC/Preproc/22_2_xcp_aroma',
        #path_output='D:/atiroms/MRI/pnTTC/pnTTC1_rsfMRI_C/31_2_conn_aroma'
        path_input='/media/veracrypt1/MRI/pnTTC/Preproc/28_2_xcp_acompcor',
        path_output='/media/veracrypt2/MRI/pnTTC/pnTTC2_rsfMRI_C/15_2_xcpout_acompcor'
        ):

        print('Starting NIfTI extraction.')

        # Create experiment folder
        print('Starting to create experiment folder.')
        list_paths_mkdir=[]
        list_paths_mkdir.append(path_output)
        list_paths_mkdir.append(os.path.join(path_output,'output'))
        list_paths_mkdir.append(os.path.join(path_output,'output','norm'))
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

        # Copy NifTI files
        print('Starting to copy NIfTI files.')
        list_dir_thread = os.listdir(os.path.join(path_input,'output'))
        list_dir_thread.sort()
        for dir_thread in list_dir_thread:
            list_dir_all=os.listdir(os.path.join(path_input,'output',dir_thread))
            list_dir_sub=[d for d in list_dir_all if os.path.isdir(os.path.join(path_input,'output',dir_thread,d)) and d.startswith('sub-')]
            list_dir_sub.sort()
            for dir_sub in list_dir_sub:
                path_file_from=dir_sub + '_img_sm6Std.nii.gz'
                path_file_from=os.path.join(path_input,'output',dir_thread,dir_sub,'norm',path_file_from)
                path_folder_to=os.path.join(path_output,'output','norm')
                shutil.copy(path_file_from,path_folder_to)
                print('Copied: '+dir_sub)
        
        print('Finished copying NIfTI files.')


##################################################
# Extract XCP-processed FC or TS data
##################################################
# This class is no longer used.
# R extract_xcp() function is used instead
class ExtractXCP():
    def __init__(self,
        path_input='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/30_xcp_36p',
        #path_output='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/42_fc',
        path_output='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/43_ts',
        atlases=['aal116','glasser360','gordon333','power264','schaefer100','schaefer200','schaefer400'],
        #suffix_result='.net',
        suffix_result='_ts.1D',
        file_id='id_5sub.txt'
        ):

        print('Starting XCP result extraction.')

        # Create output folder
        print('Starting to create output folder.')
        list_paths_mkdir=[]
        list_paths_mkdir.append(path_output)
        list_paths_mkdir.append(os.path.join(path_output,'output'))
        for p in list_paths_mkdir:
            if not os.path.exists(p):
                os.makedirs(p)
        print('Finished creating output folder.')

        # Copy log folder
        print('Starting to copy log folder.')
        path_log_in=os.path.join(path_input,'log')
        path_log_out=os.path.join(path_output,'log')
        shutil.copytree(path_log_in,path_log_out)
        print('Finished copying log folder.')

        # Copy result
        print('Starting to copy result files.')
        path_file_id=os.path.join(path_output,'log',file_id)
        with open(path_file_id, 'r') as list_id:
            list_id=list_id.readlines()
            list_id=[int(x.strip('\n')) for x in list_id]
            list_id.sort()
        for index in list_id:
            label_sub='sub-'+str(index).zfill(5)
            for atlas in atlases:
                file_from=label_sub+'_'+atlas+suffix_result
                path_file_from=os.path.join(path_input,'output',label_sub,'fcon',atlas,file_from)
                path_file_to=os.path.join(path_output,'output',file_from)
                shutil.copy(path_file_from,path_file_to)
            print('Copied results for '+ label_sub)
        print('Finished copying result files')

        print('Finishd XCP result extraction.')

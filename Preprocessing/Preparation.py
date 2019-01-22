#############
# LIBRARIES #
#############

import os
import shutil
import pandas as pd
import csv
import nilearn.image as nl_image
import json


####################################################
# INSERT SLICE TIMING AND PHASE ENCODING DIRECTION #
####################################################
# fMRIPrep preparation
# Insert slice timing and phase encoding direction data to BIDS JSON file to use in fMRIPrep

class InsertST_PED():
    def __init__(self,
        TR=2.5,
        n_slices=40,
        PED='j-',
        path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/14_bids_ses1_t1exist_boldexist/output',
        sessions=['ses-01','ses-02']
        ):

        list_dir_all = os.listdir(path_exp)
        list_dir_all.sort()
        list_dir_sub=[]
        list_dir_func=[]
        list_slicetiming=[]
        for i in range(n_slices):
            list_slicetiming.append(i*TR/n_slices)
        for dir_sub in list_dir_all:
            if dir_sub.startswith('sub-'):
                list_dir_sub.append(dir_sub)
                for session in sessions:
                    dir_func=path_exp +'/' + dir_sub + '/' + session + '/func'
                    if os.path.exists(dir_func):
                        list_dir_func.append(dir_func)
                        filename_json = dir_sub + '_' + session + '_task-rest_bold.json'
                        with open(dir_func + '/' + filename_json) as file_json_input:  
                            data = json.load(file_json_input)
                        data['SliceTiming']=list_slicetiming
                        data['PhaseEncodingDirection']=PED
                        with open(dir_func + '/' + filename_json, 'w') as file_json_output:  
                            json.dump(data, file_json_output,indent=2, sort_keys=True)
                        print('Added ST and PED data to ' + filename_json + '.')
        print('All done.')


########################
# SUBSET BIDS SUBJECTS #
########################
# fMRIPrep preparation
# Subset BIDS subjects according to available sessions or scans

class SubsetBIDS():
    def __init__(self,
        ses_remain={'ses-01','ses-02'}, # sessions not deleted 
        delete_T1only=True,
        #path_exp='/media/veracrypt1/MRI/pnTTC/BIDS/09_boldexist'
        path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/14_bids_ses1_t1exist_boldexist/output'
        ):

        list_dir_all = os.listdir(path_exp)
        list_dir_all.sort()
        for dir_sub in list_dir_all:
            if dir_sub.startswith('sub-'):
                path_sub=path_exp +'/' + dir_sub
                list_dir_ses = os.listdir(path_sub)
                for dir_ses in list_dir_ses:
                    path_ses=path_sub +'/' + dir_ses
                    delete_ses=False
                    if delete_T1only:
                        list_dir_modality=os.listdir(path_ses)
                        if not 'func' in list_dir_modality:
                            delete_ses=True
                    if not dir_ses in ses_remain:
                        delete_ses=True
                    if delete_ses:
                        shutil.rmtree(path_ses)
                    
                if len(os.listdir(path_sub))==0:
                    shutil.rmtree(path_sub)
                    print('Deleted ' + dir_sub + ' as no data exists for the subject.')

        list_dir_postremoval=os.listdir(path_exp)
        list_dir_postremoval.sort()
        list_dir_sub=[]
        for directory in list_dir_postremoval:
            if directory.startswith('sub-'):
                list_dir_sub.append(directory)

        path_tsv_original=path_exp+ '/participants.tsv'
        path_tsv_old=path_exp+ '/.participants_old.tsv'
        os.rename(path_tsv_original,path_tsv_old)

        with open(path_tsv_old,'r') as tsvin, open(path_tsv_original,'w') as tsvout:
            tsvin = csv.reader(tsvin, delimiter='\t')
            tsvout = csv.writer(tsvout, delimiter='\t')
            cnt_row=0
            for row in tsvin:
                if cnt_row==0:
                    tsvout.writerows([row])
                else:
                    if row[0] in list_dir_sub:
                        tsvout.writerows([row])
                cnt_row+=1


#######################
# SUBSET BIDS VOLUMES #
#######################
# fMRIPrep preparation
# Remove initial volumes from BIDS data

class SubsetVolume():
    def __init__(self,
        n_removevol=10,
        #path_exp='/media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/14_removeinitial'
        #path_exp='/media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/09_removeinitial'
        path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/15_bids_ses1_t1exist_boldexist_remini/output'        
        ):

        list_dir_all = os.listdir(path_exp)
        list_dir_all.sort()
        for dir_sub in list_dir_all:
            if dir_sub.startswith('sub-'):
                path_sub=path_exp +'/' + dir_sub
                list_dir_ses = os.listdir(path_sub)
                for dir_ses in list_dir_ses:
                    path_ses=path_sub +'/' + dir_ses
                    list_dir_mod=os.listdir(path_ses)
                    if 'func' in list_dir_mod:
                        file_img_in=dir_sub + '_' + dir_ses + '_task-rest_bold.nii.gz'
                        path_img=os.path.join(path_ses,'func',file_img_in)
                        img_in=nl_image.load_img(path_img)
                        img_out=img_in.slicer[:,:,:,n_removevol:]
                        img_out.to_filename(path_img)
                        print('Removed initial ' + str(n_removevol) + ' images from ' + file_img_in + '.')
        print('All done.')


###################################
# PICUP AND COPY FREESURFER FILES #
###################################
# fMRIPrep preparation
# Pickup and copy FreeSurfer-processed file to use with fMRIPrep.

class Fs2Fmriprep():
    def __init__(self,
        path_file_id='/media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/id_sub.txt',
        path_in='/media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/10_recon',
        path_out='/media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/11_fs2fmriprep'
        ):

        with open(path_file_id, 'r') as list_id:
            list_id=list_id.readlines()
            list_id=[int(x.strip('\n')) for x in list_id]
            list_id.sort()
        for i in list_id:
            path_folder_in=os.path.join(path_in,str(i).zfill(5))
            path_folder_out=os.path.join(path_out,'sub-'+str(i).zfill(5))
            shutil.copytree(path_folder_in,path_folder_out)
            print('Copied and renamed '+ path_folder_in + '.')
        print('All done.')


###################
# XCP COHORT FILE #
###################
# XCP prepataion
# Create cohort file required for xcp.

class CreateCohortfile():
    def __init__(self,
        #path_out='C:/Users/atiro/Dropbox/MRI/XCP_tutorial',
        #path_file_id='C:/Users/atiro/Dropbox/MRI/XCP_tutorial/id.txt',
        path_out='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/22_xcp_aroma_aromain/input',
        path_file_id='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/22_xcp_aroma_aromain/input/id.txt',
        #suffix_file='_ses-01_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'
        suffix_file='_ses-01_task-rest_space-T1w_desc-preproc_bold.nii.gz',
        #dir_input='10_remini_syn_12dof'
        #dir_input='16_fmriprep_newfs'
        ):

        with open(path_file_id, 'r') as list_id:
            list_id=list_id.readlines()
            list_id=[int(x.strip('\n')) for x in list_id]
            list_id.sort()
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
                                                      'input/fmriprep/sub-'+str(index).zfill(5)+'/ses-01/func/sub-'+str(index).zfill(5)+suffix_file],
                                                     index=output_func.columns),
                                           ignore_index=True)
        #output_anat.to_csv(os.path.join(path_out,'anat_cohort.csv'),index=False)
        output_func.to_csv(os.path.join(path_out,'func_cohort.csv'),index=False)
        print('All done.')


###############################
# MOVE ANAT FOLDER BEFORE XCP #
###############################
# XCP preparation
# move anat files in freesurfer output as workaround of xcp file reading error.

class MoveAnat():
    def __init__(self,
        path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/22_xcp_aroma_aromain/input/fmriprep'
        ):

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


############################
# CHANGE FOLDER PERMISSION #
############################

class FolderPermission():
    def __init__(self,
        path_input='/media/veracrypt1/MRI/pnTTC/Preproc/11_bids_ses2_t1exist',
        #path_input='/media/veracrypt1/MRI/pnTTC/Preproc/test',
        permission=0o775
        ):
        #oct(os.stat(path_input).st_mode)[-3:]
        print('File: '+ path_input)
        print('Permission: '+oct(os.stat(path_input).st_mode)[-3:])
        #print('Read permission:    '+str(os.access(path_input, os.R_OK)))
        #print('Write permission:   '+str(os.access(path_input, os.W_OK)))
        #print('Execute permission: '+str(os.access(path_input, os.X_OK)))

        for root, dirs, files in os.walk(path_input, topdown=False):
            for dir in [os.path.join(root,d) for d in dirs]:
                os.chmod(dir, permission)
            for file in [os.path.join(root, f) for f in files]:
                os.chmod(file, permission)
        
        print('All done.')

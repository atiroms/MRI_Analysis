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
# Insert information to BIDS json sidecar
##################################################
# fMRIPrep preparation
# Insert slice timing, phase encoding direction, effective echo spacing and total readout time
# data to BIDS JSON file to use in fMRIPrep

class EditJson():
    def __init__(self,
        TR=2.5,
        n_slices=40,
        PED='j-',
        EES=0.00070302532,     # Fieldmap parameter
        TRT=0.04218151959,     # Fieldmap parameter
        #path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/14_bids_ses1_t1exist_boldexist/output',
        #path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/24_st_ped/output',
        path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/test_1sub/34_bids/output',
        #sessions=['ses-01','ses-02']
        sessions=['ses-01']
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
                        data['EffectiveEchoSpacing']=EES
                        data['TotalReadoutTime']=TRT
                        with open(dir_func + '/' + filename_json, 'w') as file_json_output:  
                            json.dump(data, file_json_output,indent=2, sort_keys=True)
                        print('Modified JSON file ' + filename_json + '.')
        print('All done.')


##################################################
# Subset BIDS subjects
##################################################
# fMRIPrep preparation
# Subset BIDS subjects according to available sessions or scans

class SubsetBIDS():
    def __init__(self,
        #path_exp='/media/veracrypt1/MRI/pnTTC/BIDS/09_boldexist'
        #path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/14_bids_ses1_t1exist_boldexist/output',
        path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/19_2_fmriprep/input',
        path_file_id='/media/veracrypt1/MRI/pnTTC/Preproc/19_2_fmriprep/log/id_mild_2.csv',
        ses_remain={'ses-01','ses-02'}, # sessions not deleted 
        delete_T1only=True
        ):

        if path_file_id!="": 
            with open(path_file_id, 'r') as list_id:
                list_id=list_id.readlines()
                list_id=[int(x.strip('\n')) for x in list_id]
                list_id.sort()

            list_dir_all = os.listdir(path_exp)
            list_dir_all.sort()
            for dir_sub in list_dir_all:
                if dir_sub.startswith('sub-'):
                    id_sub=int(dir_sub.replace('sub-',''))
                    if not id_sub in list_id:
                        path_sub=path_exp +'/' + dir_sub
                        shutil.rmtree(path_sub)
                        print('Deleted ' + dir_sub + ' as subject not on the list.')
        else:
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


##################################################
# Subset BIDS volumes
##################################################
# fMRIPrep preparation
# Remove initial volumes from BIDS data

class SubsetVolume():
    def __init__(self,
        n_removevol=10,
        #path_exp='/media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/14_removeinitial'
        #path_exp='/media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/09_removeinitial'
        path_exp='/media/atiroms/MORITA_HDD4/MRI/pnTTC/Preproc/14_bids_ses1_t1exist_boldexist/output'        
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


##################################################
# Pickup and copy FreeSurfer files
##################################################
# fMRIPrep preparation
# Pickup and copy FreeSurfer-processed file to use with fMRIPrep.

class Fs2Fmriprep():
    def __init__(self,
        #path_file_id='/media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/id_sub.txt',
        #path_file_id='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/25_fmriprep/input/id_5sub.txt',
        #path_file_id='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/26_fmriprep_latest/input/id_5sub.txt',
        path_file_id='/media/veracrypt1/MRI/pnTTC/Preproc/19_2_fmriprep/log/id_mild_2.csv',
        #path_in='/media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/10_recon',
        #path_in='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/pnTTC1_T1_C_FS_10_recon/freesurfer',
        path_in='/media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/12_recon_t1exist/output',
        #path_out='/media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/11_fs2fmriprep'
        #path_out='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/25_fmriprep/output/freesurfer'
        #path_out='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/26_fmriprep_latest/output/freesurfer'
        path_out='/media/veracrypt1/MRI/pnTTC/Preproc/19_2_fmriprep/output/freesurfer'
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
        path_folder_in=os.path.join(path_in,'fsaverage')
        path_folder_out=os.path.join(path_out,'fsaverage')
        shutil.copytree(path_folder_in,path_folder_out)
        print('Copied '+ path_folder_in + '.')
        print('All done.')


##################################################
# Change folder permission
##################################################
# !DOES NOT WORK!

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

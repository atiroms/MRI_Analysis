##################################################
# Libraries
##################################################

import os
import shutil
import pandas as pd
import csv
import nilearn.image as nl_image
import nilearn.plotting as nl_plotting
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
# Calculate Phase difference Fieldmap image
##################################################
# calculate phase difference image from HeuDiConv-converted BIDS folder.
# original Philips data contains two real and two imaginary fieldmap images.

class PhaseDiff():
    def __init__(self,
        path_in='C:/Users/atiro/Dropbox/Temp/Preproc/test_1sub/32_heudiconv',
        path_out='C:/Users/atiro/Dropbox/Temp/Preproc/test_1sub/40_fieldmap',
        TE1=0.00492,
        TE2=0.00738,
        data_json_template={
            "AcquisitionMatrixPE": 64,"BodyPartExamined": "BRAIN","ConversionSoftware": "dcm2niix",
            "ConversionSoftwareVersion": "v1.0.20180622 GCC6.3.0","EchoTrainLength": 2,"FlipAngle": 60,
            "ImageOrientationPatientDICOM": [1,0,0,0,1,0],"ImageType": ["ORIGINAL","PRIMARY","P","FFE",],
            "InPlanePhaseEncodingDirectionDICOM": "COL","MRAcquisitionType": "2D",
            "MagneticFieldStrength": 3,"Manufacturer": "Philips","ManufacturersModelName": "Achieva",
            "Modality": "MR","PatientPosition": "HFS","PhaseEncodingAxis": "j",
            "PhaseEncodingSteps": 64,"ReconMatrixPE": 64,"RepetitionTime": 0.488,"SliceThickness": 4,
        }
    ):

        print('Starting PhaseDiff()')
        
        # Create output folder
        print('Starting to create output folder.')
        list_path_mkdir=[]
        list_path_mkdir.append(path_out)
        list_path_mkdir.append(os.path.join(path_out,'input'))
        for p in list_path_mkdir:
            if not os.path.exists(p):
                os.makedirs(p)
        print('Finished creating output folder.')

        # Copy log folder
        print('Starting to copy log folder.')
        path_log_in=os.path.join(path_in,'log')
        path_log_out=os.path.join(path_out,'log')
        shutil.copytree(path_log_in,path_log_out)
        print('Finished copying log folder.')

        # copy metadata files
        list_dir_subj=os.listdir(os.path.join(path_in,'output'))
        list_file_meta=[d for d in list_dir_subj if not d.startswith('sub-')]
        for f in list_file_meta:
            path_file_meta=os.path.join(path_in,'output',f)
            if os.path.isdir(path_file_meta):
                shutil.copytree(path_file_meta,os.path.join(path_out,'input',f))
            elif os.path.isfile(path_file_meta):
                shutil.copy(path_file_meta,os.path.join(path_out,'input'))
        list_dir_subj=[d for d in list_dir_subj if d.startswith('sub-')]
        list_dir_subj.sort()

        # calculate fieldmap phase difference
        print('Starting to calculate phase difference.')
        for dir_subj in list_dir_subj:
            list_dir_ses=os.listdir(os.path.join(path_in,'output',dir_subj))
            list_dir_ses.sort()
            for dir_ses in list_dir_ses:
                list_dir_seq=os.listdir(os.path.join(path_in,'output',dir_subj,dir_ses))
                list_dir_seq.sort()
                if ('anat' in list_dir_seq) and ('func' in list_dir_seq) and ('fmap' in list_dir_seq):
                    list_path_mkdir=[]
                    list_path_mkdir.append(os.path.join(path_out,'input',dir_subj))
                    list_path_mkdir.append(os.path.join(path_out,'input',dir_subj,dir_ses))
                    list_path_mkdir.append(os.path.join(path_out,'input',dir_subj,dir_ses,'fmap'))
                    for p in list_path_mkdir:
                        if not os.path.exists(p):
                            os.makedirs(p)
                    path_file_in_common=dir_subj+'_'+dir_ses+'_'
                    path_file_in_common=os.path.join(path_in,'output',dir_subj,dir_ses,
                                                    'fmap',path_file_in_common)
                    fmap_imag1=nl_image.load_img(path_file_in_common+'fieldmap1.nii.gz')
                    fmap_real1=nl_image.load_img(path_file_in_common+'fieldmap2.nii.gz')
                    fmap_imag2=nl_image.load_img(path_file_in_common+'fieldmap3.nii.gz')
                    fmap_real2=nl_image.load_img(path_file_in_common+'fieldmap4.nii.gz')

                    fmap_mag1=nl_image.math_img('np.sqrt(img1**2.0 + img2**2.0)',img1=fmap_real1,img2=fmap_imag1)
                    fmap_pha1=nl_image.math_img('np.arctan2(img2,img1)',img1=fmap_real1,img2=fmap_imag1)
                    fmap_mag2=nl_image.math_img('np.sqrt(img1**2.0 + img2**2.0)',img1=fmap_real2,img2=fmap_imag2)
                    fmap_pha2=nl_image.math_img('np.arctan2(img2,img1)',img1=fmap_real2,img2=fmap_imag2)
                    fmap_sin1=nl_image.math_img('np.sin(img)',img=fmap_pha1)
                    fmap_cos1=nl_image.math_img('np.cos(img)',img=fmap_pha1)
                    fmap_sin2=nl_image.math_img('np.sin(img)',img=fmap_pha2)
                    fmap_cos2=nl_image.math_img('np.cos(img)',img=fmap_pha2)
                    # calculate phase difference (use arctan2, sin, cos so that the range is from -pi to +pi)
                    fmap_phadiff=nl_image.math_img('np.arctan2(img2*img3-img1*img4,img2*img4+img1*img3)',
                                                   img1=fmap_sin1,img2=fmap_cos1,img3=fmap_sin2,img4=fmap_cos2)
                    # correct nifti header information (10=2+8 meaning spacial unit=mm and temporal unit=sec)
                    fmap_phadiff.header['xyzt_units']=np.array(10,dtype='uint8')
                    fmap_mag1.header['xyzt_units']=np.array(10,dtype='uint8')
                    fmap_mag2.header['xyzt_units']=np.array(10,dtype='uint8')

                    path_file_out_common=dir_subj+'_'+dir_ses+'_'
                    path_file_out_common=os.path.join(path_out,'input',dir_subj,dir_ses,'fmap',path_file_out_common)
                    fmap_phadiff.to_filename(path_file_out_common+'phasediff.nii.gz')
                    fmap_mag1.to_filename(path_file_out_common+'magnitude1.nii.gz')
                    fmap_mag2.to_filename(path_file_out_common+'magnitude2.nii.gz')

                    # Create new JSON file for phase difference image
                    data_json=data_json_template
                    data_json['EchoTime1']=TE1
                    data_json['EchoTime2']=TE2
                    data_json['IntendedFor']=dir_ses+'/func/'+dir_subj+'_'+dir_ses+'_task-rest_bold.nii.gz'
                    with open(path_file_out_common+'phasediff.json', 'w') as file_json_out:
                        json.dump(data_json, file_json_out,indent=2, sort_keys=True)

                    # copy 'anat' and 'func' folders
                    shutil.copytree(os.path.join(path_in,'output',dir_subj,dir_ses,'anat'),
                                    os.path.join(path_out,'input',dir_subj,dir_ses,'anat'))
                    shutil.copytree(os.path.join(path_in,'output',dir_subj,dir_ses,'func'),
                                    os.path.join(path_out,'input',dir_subj,dir_ses,'func'))

                    # create new '_scans.tsv' file
                    df_scans=pd.read_csv(os.path.join(path_in,'output',dir_subj,dir_ses,dir_subj+'_'+dir_ses+'_scans.tsv'),sep='\t')
                    df_scans.loc[df_scans['filename']=='fmap/'+dir_subj+'_'+dir_ses+'_fieldmap1.nii.gz','filename']='fmap/'+dir_subj+'_'+dir_ses+'_phasediff.nii.gz'
                    df_scans.loc[df_scans['filename']=='fmap/'+dir_subj+'_'+dir_ses+'_fieldmap2.nii.gz','filename']='fmap/'+dir_subj+'_'+dir_ses+'_magnitude1.nii.gz'
                    df_scans.loc[df_scans['filename']=='fmap/'+dir_subj+'_'+dir_ses+'_fieldmap3.nii.gz','filename']='fmap/'+dir_subj+'_'+dir_ses+'_magnitude2.nii.gz'
                    df_scans=df_scans.loc[df_scans['filename']!='fmap/'+dir_subj+'_'+dir_ses+'_fieldmap4.nii.gz']
                    df_scans.to_csv(os.path.join(path_out,'input',dir_subj,dir_ses,dir_subj+'_'+dir_ses+'_scans.tsv'),sep='\t',na_rep='n/a',index=False)

                    print('Calculated phase difference for subject: '+dir_subj+', ses: '+dir_ses)
                else:
                    print('Sequence missing for subj: '+dir_subj+', ses: '+dir_ses)

        # reset participants.tsv file
        list_dir_subj=os.listdir(os.path.join(path_out,'input'))
        list_dir_subj=[d for d in list_dir_subj if d.startswith('sub-')]
        list_dir_subj.sort()
        print('Number of subjects with anat, func and fmap data: '+str(len(list_dir_subj)))
        df_participants=pd.DataFrame(data={'participant_id':list_dir_subj,'age':'N/A','sex':'None','group':'control'})
        path_file_participants=os.path.join(path_out,'input','participants.tsv')
        os.remove(path_file_participants)
        df_participants.to_csv(path_file_participants,sep='\t',index=False)
        
        print('Finished PhaseDiff()')


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
        EES=0.00070302532,     # Fieldmap parameter, Brainvoyager implementation
        TRT=0.04218151959,     # Fieldmap parameter, Brainvoyager implementation
        path_exp='C:/Users/atiro/Dropbox/Temp/Preproc/test_1sub/40_fieldmap'
        ):

        print('Starting EditJson()')
        list_dir_subj=os.listdir(os.path.join(path_exp,'input'))
        list_dir_subj=[d for d in list_dir_subj if d.startswith('sub-')]
        list_dir_subj.sort()
        list_slicetiming=[i*TR/n_slices for i in range(n_slices)]
        for dir_subj in list_dir_subj:
            list_dir_ses=os.listdir(os.path.join(path_exp,'input',dir_subj))
            list_dir_ses.sort()
            for dir_ses in list_dir_ses:
                dir_func=path_exp +'/input/' + dir_subj + '/' + dir_ses + '/func'
                if os.path.exists(dir_func):
                    filename_json = dir_subj + '_' + dir_ses + '_task-rest_bold.json'
                    with open(dir_func + '/' + filename_json) as file_json_input:  
                        data = json.load(file_json_input)
                    data['SliceTiming']=list_slicetiming
                    data['PhaseEncodingDirection']=PED
                    data['EffectiveEchoSpacing']=EES
                    data['TotalReadoutTime']=TRT
                    with open(dir_func + '/' + filename_json, 'w') as file_json_output:  
                        json.dump(data, file_json_output,indent=2, sort_keys=True)
                    print('Modified JSON file ' + filename_json + '.')
        print('Finished EditJson().')


##################################################
# Subset BIDS subjects
##################################################
# fMRIPrep preparation
# Subset BIDS subjects according to available sessions or scans

class SubsetBIDS():
    def __init__(self,
        #path_exp='/media/veracrypt1/MRI/pnTTC/BIDS/09_boldexist'
        #path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/14_bids_ses1_t1exist_boldexist/output',
        path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/26_2_fmriprep/input',
        path_file_id='/media/veracrypt1/MRI/pnTTC/Preproc/26_2_fmriprep/log/w2_id_mild_2.csv',
        ses_remain={'ses-01','ses-02'}, # sessions not deleted 
        delete_T1only=True
        ):

        print('Starting SubsetBIDS()')
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

        print('Finished SubsetBIDS()')


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
        path_file_id='/media/veracrypt1/MRI/pnTTC/Preproc/26_1_fmriprep/log/w2_id_mild_1.csv',
        #path_in='/media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/10_recon',
        #path_in='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/pnTTC1_T1_C_FS_10_recon/freesurfer',
        path_in='/media/veracrypt2/MRI/pnTTC/pnTTC2_T1_C/FS/17_recon/output',
        #path_out='/media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/11_fs2fmriprep'
        #path_out='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/25_fmriprep/output/freesurfer'
        #path_out='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/26_fmriprep_latest/output/freesurfer'
        path_out='/media/veracrypt1/MRI/pnTTC/Preproc/26_1_fmriprep/output'
        ):

        print('Starting Fs2Fmriprep().')
        with open(path_file_id, 'r') as list_id:
            list_id=list_id.readlines()
            list_id=[int(x.strip('\n')) for x in list_id]
            list_id.sort()

        list_path_mkdir=[]
        list_path_mkdir.append(path_out)
        list_path_mkdir.append(os.path.join(path_out,'freesurfer'))
        for p in list_path_mkdir:
            if not os.path.exists(p):
                os.makedirs(p)
        
        for i in list_id:
            path_folder_in=os.path.join(path_in,str(i).zfill(5))
            path_folder_out=os.path.join(path_out,'freesurfer','sub-'+str(i).zfill(5))
            shutil.copytree(path_folder_in,path_folder_out)
            print('Copied and renamed '+ path_folder_in + '.')
        path_folder_in=os.path.join(path_in,'fsaverage')
        path_folder_out=os.path.join(path_out,'freesurfer','fsaverage')
        shutil.copytree(path_folder_in,path_folder_out)
        print('Copied '+ path_folder_in + '.')
        print('Finished Fs2Fmriprep().')


##################################################
# All preparation for fMRIPrep
##################################################

class PrepFmriprep():
    def __init__(self,
        #path_bids='C:/Users/atiro/Dropbox/Temp/Preproc/test_1sub/32_heudiconv',
        #path_bids='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/53_bids_fmap',
        #path_bids='/media/veracrypt1/MRI/pnTTC/Preproc/test_1sub/39_heudiconv',
        path_bids='/media/veracrypt1/MRI/pnTTC/Preproc/37_c1_bids',
        #path_bids='/media/veracrypt1/MRI/pnTTC/Preproc/38_c2_bids',
        #path_freesurfer='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/pnTTC1_T1_C_FS_10_recon',
        #path_freesurfer='/media/veracrypt1/MRI/pnTTC/Preproc/test_1sub/pnTTC1_T1_C_FS_10_recon',
        path_freesurfer='/media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/12_recon_t1exist',
        #path_freesurfer='/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/17_recon',
        #path_out='C:/Users/atiro/Dropbox/Temp/Preproc/test_1sub/40_fieldmap',
        #path_out='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/54_prep_fmriprep',
        #path_out='/media/veracrypt1/MRI/pnTTC/Preproc/test_1sub/41_prep_fmriprep',
        path_out='/media/veracrypt1/MRI/pnTTC/Preproc/39_c1_2_prep_fmriprep',
        #path_out='/media/veracrypt1/MRI/pnTTC/Preproc/40_c2_2_prep_fmriprep',
        path_file_fslicense='/usr/local/freesurfer/license.txt',
        #file_id='id_5sub.csv'
        #file_id='id_1sub.csv'
        #file_id='id_W1_T1QC_new_mild_rsfMRIexist_1.csv'
        file_id='id_W1_T1QC_new_mild_rsfMRIexist_2.csv'
        #file_id='id_W2_T1QC_new_mild_rsfMRIexist_1.csv'
        #file_id='id_W2_T1QC_new_mild_rsfMRIexist_2.csv'
        ):
        
        print('Starting PrepFmriprep()')
        _=PhaseDiff(path_in=path_bids,path_out=path_out)
        _=EditJson(path_exp=path_out)
        _=SubsetBIDS(path_exp=os.path.join(path_out,'input'),
                     path_file_id=os.path.join(path_out,'log',file_id)
                     )
        _=Fs2Fmriprep(path_in=os.path.join(path_freesurfer,'output'),
                      path_out=os.path.join(path_out,'output'),
                      path_file_id=os.path.join(path_out,'log',file_id)
                      )
        shutil.copy(path_file_fslicense,os.path.join(path_out,'log'))
        print('Finished PrepFmriprep()')


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

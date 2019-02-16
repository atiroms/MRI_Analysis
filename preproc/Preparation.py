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
# Insert slice timing and phase encoding direction
##################################################
# fMRIPrep preparation
# Insert slice timing and phase encoding direction data to BIDS JSON file to use in fMRIPrep

class InsertST_PED():
    def __init__(self,
        TR=2.5,
        n_slices=40,
        PED='j-',
        #path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/14_bids_ses1_t1exist_boldexist/output',
        path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/24_st_ped/output',
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
        suffix_file='_ses-01_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'
        #suffix_file='_ses-01_task-rest_space-T1w_desc-preproc_bold.nii.gz',
        #dir_input='10_remini_syn_12dof'
        #dir_input='16_fmriprep_newfs'
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
                                                          'input/fmriprep/sub-'+str(id_subj).zfill(5)+'/ses-01/func/sub-'+str(id_subj).zfill(5)+suffix_file],
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
        path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/33_xcp_36p_templatein/input/fmriprep'
        ):

        print('Starting to move fMRIPrep anat folder contents.')
        list_dir_all = os.listdir(path_exp)
        list_sub=[d for d in list_dir_all if os.path.isdir(os.path.join(path_exp,d)) and d.startswith('sub-')]
        list_sub.sort()
        for sub in list_sub:
            path_from=os.path.join(path_exp,sub,'anat')
            path_to=os.path.join(path_exp,sub,'ses-01','anat')
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
# joining the above two classes and some more

class XCPPrep():
    def __init__(self,
        n_proc=5,
        #path_fmriprep='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/31_fmriprep_latest_syn_templateout',
        #path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/33_xcp_36p_templatein',
        #path_fmriprep='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/31_fmriprep_latest_syn_templateout',
        #path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/37_xcp_36p_spkreg_1mm',
        #path_fmriprep='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/35_fmriprep_latest_syn_templateout_2mm',
        #path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/38_xcp_36p_spkreg_2mm',
        #path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/40_xcp_aroma_2mm',
        #path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/41_xcp_acompcor_2mm',
        #path_fmriprep='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/36_fmriprep_latest_syn_templateout_native',
        #path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/39_xcp_36p_spkreg_native',
        path_fmriprep='/media/veracrypt2/MRI/pnTTC/Preproc/test_5sub/35_fmriprep_latest_syn_templateout_2mm',
        path_exp='/media/veracrypt2/MRI/pnTTC/Preproc/test_5sub/44_xcp_parallel',
        file_id='id_5sub.txt',
        suffix_img='_ses-01_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz',
        #suffix_img='_ses-01_task-rest_space-T1w_desc-preproc_bold.nii.gz',
        #path_folder_design='/home/atiroms/Documents/GitHub/MRI_Analysis/Preprocessing/XCP_design/accessed_on_20190131/modified',
        path_folder_design='/home/atiroms/GitHub/MRI_Analysis/preproc/XCP_design/accessed_on_20190131/modified',
        #file_design='fc-36p_fconly.dsn',
        #file_design='fc-36p_fconly_old.dsn',
        #file_design='fc-36p_spkreg_fconly.dsn',
        #file_design='fc-aroma_fconly.dsn',
        #file_design='fc-acompcor_fconly.dsn',
        file_design='fc-acompcor_fconly_noqcfc.dsn',
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
                           suffix_file=suffix_img)

        # Copy fMRIPrep preprocessed data to /input folder
        print('Starting to copy fMRIPrep folder.')
        path_fmriprep_in=os.path.join(path_fmriprep,'output/fmriprep')
        path_fmriprep_out=os.path.join(path_exp,'input/fmriprep')
        shutil.copytree(path_fmriprep_in,path_fmriprep_out)
        print('Finished copying fmriprep folder.')

        # Move fmriprep /anat folder contents (workaround of XCP bug)
        _=MoveAnat(path_exp=os.path.join(path_exp,'input/fmriprep'))

        # Copy XCP design file from GitHub repsitory
        print('Starting to copy XCP design file.')
        path_dsn_in=os.path.join(path_folder_design,file_design)
        path_dsn_out=os.path.join(path_exp,'input',file_design)
        shutil.copy(path_dsn_in,path_dsn_out)
        print('Finished copying XCP design file.')

        # Generate XCP script
        print('Starting to generate XCP script.')
        _=XCPScript(n_proc, path_exp,file_design,path_img_xcp,script)
        print('Finished generating XCP script.')

        print('Finished XCP preparation.')


##################################################
# Space ID file
##################################################
# MRIQC data merging with CSUB file
# used to space skipped IDs.

class SpaceIDFile():
    def __init__(self,
        #path_file_input='/media/veracrypt1/MRI/pnTTC/Preproc/17_extractqc_ses1_t1exist/input/group_T1w.tsv',
        #path_file_output='/media/veracrypt1/MRI/pnTTC/Preproc/17_extractqc_ses1_t1exist/output/group_T1w_spaced.tsv',
        path_file_input='/media/veracrypt1/MRI/pnTTC/Preproc/17_extractqc_ses1_t1exist/input/group_bold.tsv',
        path_file_output='/media/veracrypt1/MRI/pnTTC/Preproc/17_extractqc_ses1_t1exist/output/group_bold_spaced.tsv',
        colname_id_input='bids_name',
        prefix_id_input='sub-',
        #suffix_id_input='_ses-01_T1w'
        suffix_id_input='_ses-01_task-rest_bold'
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
        df.to_csv(path_file_output,sep='\t',index=False)
        print('All done.')


##################################################
# Extract DICOM header metadata
##################################################
# used to extract scan date information

class ExtractDcmHeader():
    def __init__(self,
        path_file_in='P:/MRI/pnTTC/Preproc/00_dicom_ses12_exist/pnTTC2_T1/CSUB-00003C-02/IM-0001-0001-0001.dcm'
        ):

        img_in=pydicom.read_file(path_file_in)
        df_header=pd.DataFrame(columns=['tag','name','value'])
        for element in img_in:
            df_header=df_header.append(pd.Series([str(element.tag),element.name,element.repval],index=df_header.columns),ignore_index=True)
        self.output=df_header
        #print('DICOM element count: '+ str(len(df_header)))
        #print('All done.')

class ExtractMltDcmHeader():
    def __init__(self,
        path_exp='P:/MRI/pnTTC/Preproc/test_5sub/27_dicom_ses2_t1w/output',
        #path_exp='P:/MRI/pnTTC/Preproc/00_dicom_ses12_exist/pnTTC1_T1',
        #path_exp='P:/MRI/pnTTC/Preproc/00_dicom_ses12_exist/pnTTC2_T1',
        path_file_output='P:/MRI/pnTTC/Preproc/test_5sub/28_header_date_ses2_t1w/output/dcmheader.csv',
        #path_file_output='P:/MRI/pnTTC/Preproc/18_acquisitiondate_t1w/output/dcmheader_ses1.csv',
        #path_file_output='P:/MRI/pnTTC/Preproc/18_acquisitiondate_t1w/output/dcmheader_ses2.csv',
        keys=['Acquisition Date'],
        prefix_dir_sub='CSUB-',
        #suffix_dir_sub='C-01'
        suffix_dir_sub='C-02'
        ):

        list_dir_sub = os.listdir(path_exp)
        list_dir_sub.sort()
        df_header_mlt=pd.DataFrame(columns=['ID_pnTTC','dir_sub','key','value'])
        for dir_sub in list_dir_sub:
            ID_pnTTC=int(dir_sub.replace(prefix_dir_sub,'').replace(suffix_dir_sub,''))
            list_img=os.listdir(os.path.join(path_exp,dir_sub))
            list_img=[img for img in list_img if img.startswith('IM-')]
            list_img.sort()
            path_firstimg=os.path.join(path_exp,dir_sub,list_img[0])
            df_header_firstimg=ExtractDcmHeader(path_firstimg).output
            for key in keys:
                value=df_header_firstimg.loc[df_header_firstimg['name']==key,'value'].values[0]
                df_header_mlt=df_header_mlt.append(pd.Series([ID_pnTTC,dir_sub,key,value],index=df_header_mlt.columns),ignore_index=True)
            print('Checked directory ' + dir_sub + '.')
        
        #df_header_mlt.to_csv(path_file_output,index=False)

        list_id_input=list(df_header_mlt['ID_pnTTC'])
        list_id_output=list(range(1,max(list_id_input)+1,1))
        list_id_diff=list(set(list_id_output)-set(list_id_input))
        list_id_diff.sort()
        df_space=pd.DataFrame(np.nan,columns=df_header_mlt.columns,index=range(len(list_id_diff)))
        df_space['ID_pnTTC']=list_id_diff
        df_header_mlt=pd.concat([df_header_mlt,df_space])
        df_header_mlt=df_header_mlt.sort_values(by='ID_pnTTC')

        df_header_mlt.to_csv(path_file_output,index=False)

        print('All done.')


##################################################
# Extract XCP-processed FC or TS data
##################################################

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

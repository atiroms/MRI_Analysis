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

class PrepXCP():
    def __init__(self,
        skip_log_copy=False,
        skip_fmriprep_copy=True,
        skip_fmriprep_moveanat=True,
        n_proc=10,

        #path_fmriprep='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/401_c1_fmriprep',
        #path_exp='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/403_c1_xcp_acompcor',
        #file_id='69_id_c1_t1exist_rsfmriexist.csv',
        #ses='ses-01',
        #suffix_img='_ses-01_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz',
        #file_design='fc-acompcor_20200114.dsn',

        #path_fmriprep='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/402_c2_fmriprep',
        #path_exp='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/404_c2_xcp_acompcor',
        #file_id='68_id_c2_t1exist_rsfmriexist.csv',
        #ses='ses-02',
        #suffix_img='_ses-02_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz',
        #file_design='fc-acompcor_20200114.dsn',

        #path_fmriprep='/media/veracrypt2/MRI_img/pnTTC/preproc/401_c1_fmriprep',
        #path_exp='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/405_c1_xcp_acompcor_gsr',
        #file_id='69_id_c1_t1exist_rsfmriexist.csv',
        #ses='ses-01',
        #suffix_img='_ses-01_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz',
        #file_design='fc-acompcor_gsr_20200114.dsn',

        #path_fmriprep='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/402_c2_fmriprep',
        #path_exp='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/406_c2_xcp_acompcor_gsr',
        #file_id='68_id_c2_t1exist_rsfmriexist.csv',
        #ses='ses-02',
        #suffix_img='_ses-02_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz',
        #file_design='fc-acompcor_gsr_20200114.dsn',

        #path_fmriprep='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/401_c1_fmriprep',
        #path_exp='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/407_c1_xcp_aroma',
        #file_id='69_id_c1_t1exist_rsfmriexist.csv',
        #ses='ses-01',
        #suffix_img='_ses-01_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz',
        #file_design='fc-aroma_20200114.dsn',

        #path_fmriprep='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/402_c2_fmriprep',
        #path_exp='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/408_c2_xcp_aroma',
        #file_id='68_id_c2_t1exist_rsfmriexist.csv',
        #ses='ses-02',
        #suffix_img='_ses-02_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz',
        #file_design='fc-aroma_20200114.dsn',

        #path_fmriprep='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/401_c1_fmriprep',
        #path_exp='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/409_c1_xcp_aroma_gsr',
        #file_id='69_id_c1_t1exist_rsfmriexist.csv',
        #ses='ses-01',
        #suffix_img='_ses-01_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz',
        #file_design='fc-aroma_gsr_20200114.dsn',

        path_fmriprep='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/402_c2_fmriprep',
        path_exp='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/410_c2_xcp_aroma_gsr',
        file_id='68_id_c2_t1exist_rsfmriexist.csv',
        ses='ses-02',
        suffix_img='_ses-02_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz',
        file_design='fc-aroma_gsr_20200114.dsn',

        #path_fmriprep='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/401_c1_fmriprep',
        #path_exp='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/411_c1_xcp_36p',
        #file_id='69_id_c1_t1exist_rsfmriexist.csv',
        #ses='ses-01',
        #suffix_img='_ses-01_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz',
        #file_design='fc-36p_spkreg_20200114.dsn',

        #path_fmriprep='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/402_c2_fmriprep',
        #path_exp='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/412_c2_xcp_36p',
        #file_id='68_id_c2_t1exist_rsfmriexist.csv',
        #ses='ses-02',
        #suffix_img='_ses-02_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz',
        #file_design='fc-36p_spkreg_20200114.dsn',

        #path_folder_design='/home/atiroms/Documents/GitHub/MRI_Analysis/Preprocessing/XCP_design/accessed_on_20190131/modified',
        #path_folder_design='/home/atiroms/GitHub/MRI_Analysis/preproc/XCP_design/accessed_on_20190131/modified',
        #path_folder_design='C:/Users/NICT_WS/GitHub/MRI_Analysis/preproc/XCP_design/accessed_on_20190131/modified',
        path_folder_design='/home/atiroms/GitHub/MRI_Analysis/preproc/XCP_design/20200114_censoring',
        #path_folder_design='/home/atiroms/GitHub/MRI_Analysis/preproc/XCP_design/20200326_nocensoring',

        #path_img_xcp='/data/applications/xcpEngine-070-20190130.simg',
        #path_img_xcp='/data/applications/xcpEngine-070-20190311.simg',
        #path_img_xcp='/data/applications/xcpEngine-100-20190628.simg',
        path_img_xcp='/data/applications/xcpEngine-100-20200113.simg',
        script='export BRAINSPACE="${HOME}/data/input/space"\nexport BRAINATLAS="${HOME}/data/input/atlas"\nsingularity run -B {path_exp}:${HOME}/data {path_img_xcp} -d ${HOME}/data/input/{file_design} -c ${HOME}/data/input/func_cohort_{id_proc}.csv -o ${HOME}/data/output/{id_proc} -t 1 -r ${HOME}/data'
        #script='singularity run --cleanenv -B {path_exp}:${HOME}/data {path_img_xcp} -d ${HOME}/data/input/{file_design} -c ${HOME}/data/input/func_cohort_{id_proc}.csv -o ${HOME}/data/output/{id_proc} -t 1 -r ${HOME}/data'
        ):

        print('Starting PrepXCP().')

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
        if not skip_log_copy:
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

        print('Finished PrepXCP().')


##################################################
# Multiple XCP preparation
##################################################
# batch run of above

class MultiPrepXCP():
    def __init__(self,
        prefix_path_fmriprep='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/55_',
        suffix_path_fmriprep='_fmriprep',
        prefix_path_exp='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/56_',
        suffix_path_exp='_prestats',
        #list_iteration=['01','02','03','04','05','06']
        #list_iteration=['07','08','09','10','11','12'],
        list_iteration=['07','08','09','11'],
        ):

        print('Starting MultiPrepXCP().')
        for itr in list_iteration:
            path_fmriprep=prefix_path_fmriprep+itr+suffix_path_fmriprep
            path_exp=prefix_path_exp+itr+suffix_path_exp
            _=PrepXCP(path_fmriprep=path_fmriprep,path_exp=path_exp)
            print('Finished preparation of: '+itr)

        print('Finished MultiPrepXCP().')


##################################################
# Extract Normalized NIfTI data
##################################################
# Extract XCP-preprocessed NIfTI data

class ExtractNifti():
    def __init__(self,
        #path_input='/media/veracrypt2/MRI_img/pnTTC/preproc/371_c1_xcp_acompcor_gsr',
        #path_output='/media/veracrypt2/MRI_img/pnTTC/preproc/378_nii_acompcor_gsr',
        #ses='ses-01'

        #path_input='/media/veracrypt2/MRI_img/pnTTC/preproc/372_c2_xcp_acompcor_gsr',
        #path_output='/media/veracrypt2/MRI_img/pnTTC/preproc/378_nii_acompcor_gsr',
        #ses='ses-02'

        path_input='/media/veracrypt2/MRI_img/pnTTC/preproc/381_c1_xcp_aroma_gsr',
        path_output='/media/veracrypt2/MRI_img/pnTTC/preproc/388_nii_aroma_gsr',
        ses='ses-01'

        #path_input='/media/veracrypt2/MRI_img/pnTTC/preproc/382_c2_xcp_aroma_gsr',
        #path_output='/media/veracrypt2/MRI_img/pnTTC/preproc/388_nii_aroma_gsr',
        #ses='ses-02'

        #path_input='/media/veracrypt2/MRI_img/pnTTC/preproc/191_c1_xcp_36p',
        #path_output='/media/veracrypt2/MRI_img/pnTTC/preproc/198_nii_36p',
        #ses='ses-01'

        #path_input='/media/veracrypt2/MRI_img/pnTTC/preproc/192_c2_xcp_36p',
        #path_output='/media/veracrypt2/MRI_img/pnTTC/preproc/198_nii_36p',
        #ses='ses-02'
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
                path_file_to=ses+'_'+dir_sub+'.nii.gz'
                path_file_to=os.path.join(path_output,'output','norm',path_file_to)
                if os.path.exists(path_file_from):
                    shutil.copy(path_file_from,path_file_to)
                    print('Copied: '+dir_sub)
        
        print('Finished copying NIfTI files.')


##################################################
# Extract n*_quality.csv data
##################################################
# Extract from XCP result folder

class ExtractQuality():
    def __init__(self,
        #path_input='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/403_c1_xcp_acompcor',
        #path_output='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/423_c1_quality_acompcor',
        #path_input='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/404_c2_xcp_acompcor',
        #path_output='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/424_c2_quality_acompcor',
        #path_input='J:/MRI_img/pnTTC/preproc/405_c1_xcp_acompcor_gsr',
        #path_output='C:/Users/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP/395_c1_quality',
        #path_input='J:/MRI_img/pnTTC/preproc/406_c2_xcp_acompcor_gsr',
        #path_output='C:/Users/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP/396_c2_quality',
        path_input='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/434_c2_xcp_acompcor',
        path_output='/media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/454_c2_quality_acompcor',

        skip_mkdir=False,
        skip_copylog=False,
        filename_output='quality.csv'
        ):

        print('Starting ExtractQuality().')

        if not skip_mkdir:
            # Create experiment folder
            print('Starting to create experiment folder.')
            list_paths_mkdir=[]
            list_paths_mkdir.append(path_output)
            list_paths_mkdir.append(os.path.join(path_output,'output'))
            for p in list_paths_mkdir:
                if not os.path.exists(p):
                    os.makedirs(p)
            print('Finished creating experiment folder.')

        if not skip_copylog:
            # Copy log file
            print('Starting to copy log folder.')
            path_log_in=os.path.join(path_input,'log')
            path_log_out=os.path.join(path_output,'log')
            shutil.copytree(path_log_in,path_log_out)
            print('Finished copying log folder.')

        # read quality data
        list_dir_thread = os.listdir(os.path.join(path_input,'output'))
        list_dir_thread.sort()
        df_quality=pd.DataFrame()
        for i in range(len(list_dir_thread)):
            path_quality=os.path.join(path_input,'output',list_dir_thread[i],'group')
            path_file_quality=glob.glob(path_quality+'/n*_quality.csv')[0]
            df_quality_add=pd.read_csv(path_file_quality)
            df_quality=pd.concat([df_quality,df_quality_add])
            print('Finished extracting from thread: ',list_dir_thread[i])
        df_quality.loc[:,'id0']=[int(i.replace('sub-','')) for i in df_quality.loc[:,'id0']]
        df_quality_spaced=pd.DataFrame([i for i in range(1,max(df_quality.loc[:,'id0'])+1)],columns=['id0'])
        df_quality_spaced=pd.merge(df_quality_spaced,df_quality,how='left',on='id0')
        df_quality_spaced=df_quality_spaced.rename(columns={'id0':'ID_pnTTC'})
        path_file_output=os.path.join(path_output,'output',filename_output)
        df_quality_spaced.to_csv(path_file_output,index=False)

        print('Finished ExtractQuality().')


##################################################
# Multiple Extract n*_quality.csv data
##################################################
# batch run of above

class MultiExtractQuality():
    def __init__(self,
        prefix_path_input='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/56_',
        suffix_path_input='_prestats',
        path_output='/media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/57_quality',
        #list_iteration=['01','02','03','04','05','06'],
        #list_iteration=['07','08','09','10','11','12'],
        list_iteration=['07','08','09','11'],
        ):

        print('Starting MultiExtractQuality().')

        print('Starting to create experiment folder.')
        list_paths_mkdir=[]
        list_paths_mkdir.append(path_output)
        list_paths_mkdir.append(os.path.join(path_output,'output'))
        list_paths_mkdir.append(os.path.join(path_output,'output','quality'))
        for p in list_paths_mkdir:
            if not os.path.exists(p):
                os.makedirs(p)
        print('Finished creating experiment folder.')

        # Copy log file (only for the first input folder)
        print('Starting to copylog folder.')
        path_log_in=os.path.join(prefix_path_input+list_iteration[0]+suffix_path_input,'log')
        path_log_out=os.path.join(path_output,'log')
        shutil.copytree(path_log_in,path_log_out)
        print('Finished copying log folder.')

        for itr in list_iteration:
            path_input=prefix_path_input+itr+suffix_path_input
            filename_output=itr+'_quality.csv'
            _=ExtractQuality(path_input=path_input,path_output=path_output,
                             filename_output=filename_output,
                             skip_mkdir=True,skip_copylog=True)
            print('Finished extracting quality for '+itr)
        
        print('Finished MultiExtractQuality().')


##################################################
# Extraction of NIfTI and quality data
##################################################

class PostXCP():
    def __init__(self,
        path_input='/media/veracrypt1/MRI/pnTTC/Preproc/43_c1_2_xcp_acompcor',
        path_output='/media/veracrypt1/MRI/pnTTC/Preproc/45_c1_2_nii_acompcor'
        ):

        print('Starting PostXCP().')
        _=ExtractNifti(path_input=path_input,path_output=path_output)
        _=ExtractQuality(path_input=path_input,path_output=path_output,
                         skip_mkdir=False,skip_copylog=True)
        print('Finished PostXCP().')
        

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

##################################################
# Libraries
##################################################

import os
import math
import numpy as np
import pandas as pd
import shutil


##################################################
# Check FreeSurfer files
##################################################
# Check FreeSurfer folder for subject list, error log and folder size.

class CheckFreesurfer():
    def __init__(self,
        #path_exp='/media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/10_recon',
        #path_exp='/media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/10.1_recon_t1qcout/output',
        #path_exp='/media/atiroms/MORITA_HDD4/MRI/pnTTC/pnTTC1_T1_C/FS/12_recon_t1exist/output',
        path_exp='/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/17_recon/output',
        string_log_ok='finished without error at',
        #file_output='/media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/check_freesurfer.csv'
        file_output='/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/17_recon/log/18_checkfreesurfer.csv'
        ):

        list_dir_all = os.listdir(path_exp)
        list_sub = [d for d in list_dir_all if (os.path.isdir(os.path.join(path_exp,d)) and not d.startswith('fsaverage'))]
        list_sub.sort()
        df_out=pd.DataFrame(columns=['sub','log','size'])
        list_log_error=[]
        for sub in list_sub:
            print('Checking subject '+ str(sub)+ ' ...')
            with open(os.path.join(path_exp,sub,'scripts/recon-all.log'),'r') as file_log:
                text_log=file_log.read()
            if string_log_ok in text_log:
                log_ok=1
            else:
                log_ok=0
                list_log_error.append(sub)
            df_out=df_out.append(pd.Series([sub,log_ok,self.get_dir_size(path=os.path.join(path_exp,sub))],index=df_out.columns),ignore_index=True)
        df_out.to_csv(file_output,index=False)
        print('Total FreeSurfer subject folders: ' + str(len(list_sub)))
        print('FreeSurfer Log Error in:')
        print(list_log_error)
        print('All done.')

    def get_dir_size(self, path='.'):
        total = 0
        with os.scandir(path) as it:
            for entry in it:
                if entry.is_file():
                    total += entry.stat().st_size
                elif entry.is_dir():
                    total += self.get_dir_size(entry.path)
        return total


##################################################
# Zero-pad and concatenate id file
##################################################
# read id file, zero-pad and concatenate into a string
# for general use

class ZeropadConcat():
    def __init__(self,
        path_file_id='/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/18_meas/log/list_id.txt',
        n_zfill=5,
        path_file_output='/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/18_meas/log/str_id.txt'
        ):
        with open(path_file_id, 'r') as list_id:
            list_id=list_id.readlines()
            list_id=[int(x.strip('\n')) for x in list_id]
            list_id.sort()
        text_id=''
        for index in list_id:
            text_id=text_id + ' ' + str(index).zfill(n_zfill)
        
        with open(path_file_output,'w') as file_output:
            file_output.write(text_id)
        print('All done.')


##################################################
# read ID file into list
##################################################

class ReadID():
    def __init__(self,
        #path_file_id='/media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/script/id_11_recon.txt'
        path_file_id='/media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/10.1_recon_t1qcout/input/id.txt'
        ):

        file=open(path_file_id, 'r')
        file=file.readlines()
        self.output=[int(x.strip('\n')) for x in file]


##################################################
# Scan FreeSurfer folder into ID list
##################################################

class ScanFSFolder():
    def __init__(self,
        #path_exp='/media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/10_recon',
        #path_exp='/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/15_recon',
        path_exp='/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/17_recon/output',
        list_exceptions=['fsaverage', 'id.txt','script.txt']
        ):

        list_dir = os.listdir(path_exp)
        list_dir.sort()
        #self.output = [int(i) for i in list_dir if i != 'fsaverage' and i != 'id.txt' and i != 'script.txt']
        self.output = [int(i) for i in list_dir if i not in list_exceptions]


##################################################
# Scan .nii.gz folder into ID list
##################################################

class ScanNiiFolder():
    def __init__(self,
        path_exp='/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/13_nii.gz_unite',
        list_exceptions=['fsaverage', 'id.txt','script.txt']
        ):

        list_dir = os.listdir(path_exp)
        list_dir.sort()
        list_id= [int(i.replace('CSUB-','').replace('C-02.nii.gz','')) for i in list_dir if i not in list_exceptions]
        list_id.sort()
        self.output=list_id


##################################################
# save list id into text file
##################################################

class SaveListID():
    def __init__(self,
        list_id,
        #path_file_output='/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/16_recon_t1qcout/log/list_id.txt'
        path_file_output='/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/18_meas/log/list_id.txt'
        ):

        self.output=''
        for item in list_id:
            self.output=self.output + str(item) + '\n'
        file=open(path_file_output,'w')
        file.write(self.output)
        file.close()


##################################################
# generate script from ID file
##################################################

class GenerateScript():
    def __init__(self,
        list_id,
        path_src='',
        path_dst='',
        id_ses=0,
        head_script='SUBJECTS_DIR={path_src}/output\ncd $SUBJECTS_DIR\n',
        script='recon-all -i {path_dst}/output/sub-{id_sub}_ses-{id_ses}_T1w.nii -subject {id_sub} -all -qcache',
        connector=' ; '
        ):
        self.output=head_script
        for id_subj in list_id:
            id_sub=str(id_subj).zfill(5)
            id_ses=str(id_ses).zfill(2)
            script_subj=script
            script_subj=script_subj.replace('{path_src}',path_dst).replace('{path_dst}',path_dst)
            script_subj=script_subj.replace('{id_sub}',id_sub).replace('{id_ses}',id_ses)
            self.output=self.output+script_subj
            if id_subj != list_id[-1]:
                self.output=self.output + connector


##################################################
# Generate FreeSurfer scripts for multiprocessing
##################################################

class PrepFS():
    def __init__(self,
        file_id='id_fs.csv',
        path_src='',
        path_dst='',
        id_ses=1,
        file_script='freesurfer_mp.sh',
        n_proc=28
        ):

        print('Starting PrepFS()')

        # Create output folder
        print('Starting to create output folder.')
        list_path_mkdir=[]
        list_path_mkdir.append(path_dst)
        list_path_mkdir.append(os.path.join(path_dst,'output'))
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

        # Create ID list
        with open(os.path.join(path_dst,'log',file_id), 'r') as list_id:
            list_id=list_id.readlines()
            list_id=[int(x.strip('\n')) for x in list_id]
            list_id.sort()
        print('Number of subjects: '+str(len(list_id)))

        # Create list of ID lists
        n_subj=len(list_id)
        n_subj_per_proc=int(np.ceil(n_subj/n_proc))
        n_proc_floor=n_proc*n_subj_per_proc-n_subj
        n_proc_ceil=n_proc-n_proc_floor
        print('Multiprocessing: '+str(n_subj)+' total subs, '+str(n_proc)+' total procs: '+str(n_proc_ceil)+' procs x '
              +str(n_subj_per_proc)+' subs, '+str(n_proc_floor)+' procs x '+ str(n_subj_per_proc-1)+' subs.')
        list_list_id=[]
        for id_proc in range(n_proc):
            if id_proc<n_proc_ceil:
                list_list_id.append(list_id[(id_proc*n_subj_per_proc):((id_proc+1)*n_subj_per_proc)])
            else:
                list_list_id.append(list_id[(n_subj_per_proc*n_proc_ceil+(id_proc-n_proc_ceil)*(n_subj_per_proc-1)):
                                            (n_subj_per_proc*n_proc_ceil+(id_proc-n_proc_ceil+1)*(n_subj_per_proc-1))])
        
        # Create multiple scripts
        script_out=''
        for id_proc in range(n_proc):
            list_id=list_list_id[id_proc]
            str_id_proc=str(id_proc+1).zfill(2)
            script_proc_head='# Process No. '+ str_id_proc + '\n'
            script_proc=GenerateScript(list_id=list_id,path_src=path_src,path_dst=path_dst,id_ses=id_ses).output
            script_proc_tail='\n\n'
            script_proc=script_proc_head+script_proc+script_proc_tail
            script_out=script_out+script_proc

        file=open(os.path.join(path_dst,'log'),'w')
        file.write(script_out)
        file.close()


##################################################
# generate multiple scripts for remaining subjects
##################################################
# !No longer used!

#class RemainingScript():
#    def __init__(self):
#        list_diff=list(set(ReadID().output)-set(ScanFSFolder().output))
#        GenerateMultiScript(list_diff)
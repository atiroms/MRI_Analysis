##################################################
# Libraries
##################################################

import os
import math


##################################################
# Zero-pad and concatenate id file
##################################################
# read id file, zero-pad and concatenate into a string
# for general use

class ZeropadConcat():
    def __init__(self,
        path_file_id='/media/veracrypt1/MRI/pnTTC/BIDS/misc/id_W2_T1exist.txt',
        n_zfill=5,
        path_file_output='/media/veracrypt1/MRI/pnTTC/BIDS/misc/id_string_W2_T1exist.txt'
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
        path_exp='/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/15_recon',
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
        path_file_output='/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/16_recon_t1qcout/log/list_id.txt'
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
        #head='SUBJECTS_DIR=/media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/10.1_recon_t1qcout/output\ncd $SUBJECTS_DIR\n',
        head='SUBJECTS_DIR=/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/16_recon_t1qcout/output\ncd $SUBJECTS_DIR\n',
        #head='SUBJECTS_DIR=/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/15_recon\ncd $SUBJECTS_DIR\n',
        #text=['recon-all -i /media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/06_qc/CSUB-',
        #      'C-01.nii -subject ',
        #      ' -all -qcache'],
        #text=['recon-all -i /media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/14_qc/CSUB-',
        #      'C-02.nii.gz -subject ',
        #      ' -all -qcache'],
        #text=['recon-all -i /media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/04_nii/CSUB-',
        #      'C-01.nii -subject ',
        #      ' -all -qcache'],
        text=['recon-all -i /media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/13_nii.gz_unite/CSUB-',
              'C-02.nii -subject ',
              ' -all -qcache'],
        connector=' ; '
        ):

        self.output=head
        for i in list_id:
            script=text[0] + str(i).zfill(5) + text[1] + str(i).zfill(5) + text[2]
            self.output=self.output + script
            if i != list_id[-1]:
                self.output=self.output + connector


##################################################
# generate scripts for manual multiprocessing
##################################################

class GenerateMultiScript():
    def __init__(self,
        list_id,
        #path_file_out='/media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/script',
        #path_file_out='/media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/10.1_recon_t1qcout/log/10.1_recon_t1qcout.sh',
        path_file_out='/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/16_recon_t1qcout/log/16_recon_t1qcout.sh',
        n_scripts=19
        ):

        len_list=len(list_id)
        len_script=math.ceil(len_list/n_scripts)
        len_remaining=len_list
        output=''
        for i in range(n_scripts):
            list_script=list_id[(i*len_script):((i+1)*len_script)]
            output=output + GenerateScript(list_id=list_script).output + '\n\n'
            len_remaining=len_list-len_script
            if len_remaining<1:
                break
        file=open(path_file_out,'w')
        file.write(output)
        file.close()


##################################################
# generate multiple scripts for remaining subjects
##################################################

class RemainingScript():
    def __init__(self):
        list_diff=list(set(ReadID().output)-set(ScanFSFolder().output))
        GenerateMultiScript(list_diff)
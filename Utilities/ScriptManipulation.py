
##############
# PARAMETERS #
##############

#path_exp='D:/MRI/pnTTC2_T1_C/FS/10_recon'
#path_exp='/media/veracrypt1/MRI/pnTTC2_T1_C/FS/12_nii.gz'
#path_exp='/media/veracrypt1/MRI/pnTTC2_T1_C/FS/11_dcm'
#path_exp='/media/veracrypt1/MRI/pnTTC2_T1_C/FS/10_tar.gz_new'
#path_exp='/media/veracrypt1/MRI/pnTTC2_rsfMRI_C/10_tar.gz_new'
#path_exp='/media/veracrypt1/MRI/pnTTC2_rsfMRI_C/11_dcm'
#path_exp='/media/veracrypt1/MRI/pnTTC2_rsfMRI_C/12_nii.gz'
path_exp = '/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/14_qc'

file_id='id.txt'


#head='SUBJECTS_DIR=/media/veracrypt1/MRI/pnTTC1_T1_C/FS/10_recon\ncd $SUBJECTS_DIR\n'
head='SUBJECTS_DIR=/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/15_recon\ncd $SUBJECTS_DIR\n'

#text=['recon-all -i /media/veracrypt1/MRI/pnTTC1_T1_C/FS/06_qc/CSUB-',
#      'C-01.nii -subject ',
#      ' -all -qcache']

text=['recon-all -i /media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/14_qc/CSUB-',
      'C-02.nii.gz -subject ',
      ' -all -qcache']

connector=' ; '

n_scripts=30


#############
# LIBRARIES #
#############

import os
import math


class Read_ID():
    def __init__(self,path_exp=path_exp, file_id=file_id):
        file=open(path_exp + '/' + file_id, 'r')
        file=file.readlines()
        self.output=[int(x.strip('\n')) for x in file]


class Extract_FolderID():
    def __init__(self,path_exp=path_exp):
        self.path_exp=path_exp
        self.list_dir = os.listdir(self.path_exp)
        #self.output = [int(i) for i in self.list_dir if i != 'fsaverage' and i != 'id.txt' and i != 'script.txt']
        self.output = [i for i in self.list_dir if i != 'fsaverage' and i != 'id.txt' and i != 'script.txt']

class Save_List_ID():
    def __init__(self,list_id,path_exp=path_exp):
        self.path_exp=path_exp
        self.list_id=list_id
        self.output=''
        for item in self.list_id:
            self.output=self.output + item + '\n'
        file=open(path_exp + '/list_id.txt','w')
        file.write(self.output)
        file.close()


class Get_Diff():
    def __init__(self):
        self.output=list(set(Read_ID().output)-set(Extract_FolderID().output))


class Generate_Script():
    def __init__(self, list_id,head=head,text=text,connector=connector):
        self.list_id=list_id
        self.head=head
        self.text=text
        self.connector=connector
        self.output=head
        for i in self.list_id:
            script=text[0] + str(i).zfill(5) + text[1] + str(i).zfill(5) + text[2]
            self.output=self.output + script
            if i != self.list_id[-1]:
                self.output=self.output + connector

class Generate_MultiScript():
    def __init__(self, list_id, n_scripts=n_scripts, path_exp=path_exp):
        self.list_id=list_id
        len_list=len(self.list_id)
        len_script=math.ceil(len_list/n_scripts)
        len_remaining=len_list
        self.output=''
        for i in range(n_scripts):
            list_script=self.list_id[(i*len_script):((i+1)*len_script)]
            self.output=self.output + Generate_Script(list_id=list_script).output + '\n\n'
            len_remaining=len_list-len_script
            if len_remaining<1:
                break
        file=open(path_exp + '/script.txt','w')
        file.write(self.output)
        file.close()


class Remaining_Script():
    def __init__(self):
        list_id=Get_Diff().output
        Generate_MultiScript(list_id)
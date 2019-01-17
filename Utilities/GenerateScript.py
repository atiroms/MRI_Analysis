#############
# LIBRARIES #
#############

import os
import math


############################
# ZERO-PAD AND CONCATENATE #
############################
# read id file, zero-pad and concatenate into a string
# for general use

class ZeropadConcat():
    def __init__(self):
        ############
        # Parameters
        path_file_id='/media/veracrypt1/MRI/pnTTC/BIDS/misc/id_W1_T1exist.txt'
        n_zfill=5
        path_file_output='/media/veracrypt1/MRI/pnTTC/BIDS/misc/id_string_W1_T1exist.txt'
        ############

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


##########################
# READ ID FILE INTO LIST #
##########################

class ReadID():
    def __init__(self):
        ############
        # Parameters
        path_file_id=''
        ############

        file=open(path_file_id, 'r')
        file=file.readlines()
        self.output=[int(x.strip('\n')) for x in file]


#########################
# SCAN FOLDER INTO LIST #
#########################

class ExtractFolderID():
    def __init__(self):
        ############
        # Parameters
        path_exp=''
        list_exceptions=['fsaverage', 'id.txt','script.txt']
        ############

        list_dir = os.listdir(path_exp)
        #self.output = [int(i) for i in list_dir if i != 'fsaverage' and i != 'id.txt' and i != 'script.txt']
        self.output = [int(i) for i in list_dir if i not in list_exceptions]


###############################
# SAVE LIST ID INTO TEXT FILE #
###############################

class SaveListID():
    def __init__(self,list_id):
        ############
        # Parameters
        path_file_output=''
        ############

        self.output=''
        for item in list_id:
            self.output=self.output + item + '\n'
        file=open(path_file_output,'w')
        file.write(self.output)
        file.close()


################################
# GENERATE SCRIPT FROM ID FILE #
################################

class GenerateScript():
    def __init__(self,list_id):
        ############
        # Parameters
        head='SUBJECTS_DIR=/media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/10_recon\ncd $SUBJECTS_DIR\n'
        #head='SUBJECTS_DIR=/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/15_recon\ncd $SUBJECTS_DIR\n'
        text=['recon-all -i /media/veracrypt1/MRI/pnTTC/pnTTC1_T1_C/FS/06_qc/CSUB-',
              'C-01.nii -subject ',
              ' -all -qcache']
        #text=['recon-all -i /media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/14_qc/CSUB-',
        #      'C-02.nii.gz -subject ',
        #      ' -all -qcache']
        connector=' ; '
        ############

        self.output=head
        for i in list_id:
            script=text[0] + str(i).zfill(5) + text[1] + str(i).zfill(5) + text[2]
            self.output=self.output + script
            if i != list_id[-1]:
                self.output=self.output + connector


###############################################
# GENERATE SCRIPTS FOR MANUAL MULTIPROCESSING #
###############################################

class GenerateMultiScript():
    def __init__(self, list_id):
        ############
        # Parameters
        path_exp=''
        #n_scripts=30
        n_scripts=60
        ############
        
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
        file=open(path_exp + '/script.txt','w')
        file.write(output)
        file.close()


####################################################
# GENERATE MULTIPLE SCRIPTS FOR REMAINING SUBJECTS #
####################################################

class RemainingScript():
    def __init__(self):
        list_diff=list(set(ReadID().output)-set(ExtractFolderID().output))
        GenerateMultiScript(list_diff)
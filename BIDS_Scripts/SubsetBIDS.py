

##############
# Parameters #
##############

subset_ses='ses-01'
subset_T1only=True
#path_exp='C:/Users/atiro/Documents/MRI/pnTTC/BIDS/test_3sub/05_subset'
path_exp='/media/veracrypt1/MRI/pnTTC/BIDS/03_ses1'

#############
# LIBRARIES #
#############

import os
import shutil
import csv
#import json


#############
# MAIN CODE #
#############

class SubsetBIDS():
    def __init__(self,path_exp=path_exp,subset_ses=subset_ses,subset_T1only=subset_T1only):
        self.path_exp=path_exp
        list_dir_all = os.listdir(self.path_exp)
        for dir_sub in list_dir_all:
            if dir_sub.startswith('sub-'):
                path_sub=self.path_exp +'/' + dir_sub
                list_dir_ses = os.listdir(path_sub)
                for dir_ses in list_dir_ses:
                    path_ses=path_sub +'/' + dir_ses
                    delete_sub=False
                    if subset_T1only:
                        list_dir_modality=os.listdir(path_ses)
                        if not 'func' in list_dir_modality:
                            delete_sub=True
                    if not dir_ses==subset_ses:
                        delete_sub=True
                    if delete_sub:
                        shutil.rmtree(path_ses)
                    
                if len(os.listdir(path_sub))==0:
                    shutil.rmtree(path_sub)
                    print('Deleted ' + dir_sub + ' as no data exists for the subject.')

        list_dir_postremoval=os.listdir(self.path_exp)
        list_dir_sub=[]
        for directory in list_dir_postremoval:
            if directory.startswith('sub-'):
                list_dir_sub.append(directory)

        path_tsv_original=self.path_exp+ '/participants.tsv'
        path_tsv_old=self.path_exp+ '/.participants_old.tsv'
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
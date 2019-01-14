#############
# LIBRARIES #
#############

import os
import shutil
import csv
#import json
import nilearn.image as nl_image


# for subsetting participants/sessions from BIDS data
# used for extracting pnTTC W1 data from W1 + W2 data 
##############
# Parameters #
##############

subset_ses={'ses-01','ses-02'}
subset_T1only=True
path_exp='/media/veracrypt1/MRI/pnTTC/BIDS/09_boldexist'

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
                    delete_ses=False
                    if subset_T1only:
                        list_dir_modality=os.listdir(path_ses)
                        if not 'func' in list_dir_modality:
                            delete_ses=True
                    if not dir_ses in subset_ses:
                        delete_ses=True
                    if delete_ses:
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



# for removing initial volumes from BIDS format data
##############
# Parameters #
##############
#n_removevol=10
#path_exp='/media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/14_removeinitial'
#path_exp='/media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/09_removeinitial'
#path_exp='/media/veracrypt1/MRI/pnTTC/BIDS/07_removeinitial'

class SubsetVolume():
    def __init__(self, path_exp=path_exp,n_removevol=n_removevol):
        list_dir_all = os.listdir(path_exp)
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

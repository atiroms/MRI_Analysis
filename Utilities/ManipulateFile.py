#############
# LIBRARIES #
#############

import os
import shutil
import glob
import numpy as numpy
import pandas as pd

###############
# PICKUP FILE #
###############
# for general file selection according to list of IDs

class PickupFile():
    def __init__(self):
        ############
        # Parameters
        path_from=''
        path_to=''
        path_id=''
        file_id=''
        prefix=''
        suffix=''
        ############

        file_id_open=open(path_id + '/' + file_id, 'r')
        file_id_open=file_id_open.readlines()
        id_list=[int(x.strip('\n')) for x in file_id_open]
        for id in id_list:
            name_file=prefix + str(id).zfill(5) + suffix
            file_from=path_from + '/' + name_file
            file_to=path_to + '/' + name_file
            shutil.copy(file_from, file_to)
            print('Copied file: ' + file_from)
        print('Finished file pick-up.')


############################
# COPY FILE NAME TO FOLDER #
############################
# copy file name to folder name and put file into folder
# used as preparation of DPARSF

class File2Folder():
    def __init__(self):
        ############
        # Parameters
        suffix = '.nii'
        #path='J:/MRI/pnTTC/Prosociality_DC_Dr_Okada/Analysis/FunImg'
        path='J:/MRI/pnTTC/Prosociality_DC_Dr_Okada/Analysis/T1Img'
        ############

        list_dir = os.listdir(path)
        for f in list_dir:
            if f.endswith(suffix):
                name_folder = f[:-(len(suffix))]
                path_newfolder = path + '/' + name_folder + '/'
                path_oldfile = path + '/' + f
                if not os.path.exists(path_newfolder):
                    os.makedirs(path_newfolder)
                shutil.move(path_oldfile, path_newfolder)
                print('Moved ' + f + '.')


###################
# PICKUP ROI FILE #
###################
#  used for sorting DPARSF-analyzed data into folder

class PickupROIFile():
    def __init__(self):
        ############
        # Parameters
        path_from='H:/MRI/pnTTC/Prosociality_DC_Dr_Okada/SPM/FC_FunRawARCWSF'
        path_to='H:/MRI/pnTTC/Prosociality_DC_Dr_Okada/SPM'
        path_list_id='H:/MRI/pnTTC/Prosociality_DC_Dr_Okada/SPM/Script'

        #list_rois=[1,2,3,4,5,6,7,8]
        list_rois=[6]
        #list_models=[1,2,3]
        list_models=[2]
        ############

        for model in list_models:
            for roi in list_rois:
                path_folder_to = path_to + '/Model' + str(model) + '/ROI' + str(roi) + '/FC'
                if not os.path.exists(path_folder_to):
                    os.makedirs(path_folder_to)
                with open(path_list_id+ '/id_model' + str(model) + '.txt', 'r') as list_id:
                    list_id=list_id.readlines()
                    list_id=[x.strip('\n') for x in list_id]
                for i in list_id:
                    path_file_from=path_from + '/zROI' + str(roi) + 'FCMap_CSUB-' + str(i).zfill(5) + 'C-01.nii'
                    if os.path.exists(path_file_from):
                        shutil.copy(path_file_from,path_folder_to)
                    else:
                        print('File ' + path_file_from + ' does not exist.')
                print('Copied Model ' + str(model) + ', ROI ' + str(roi) + ' files.')
        print('All done.')


#################
# SORT ROI FILE #
#################
# used for sorting DPARSF-analyzed ROI FC file into groups

class SortROIFile():
    def __init__(self):
        ############
        # Parameters
        path_from='H:/MRI/pnTTC/Prosociality_DC_Dr_Okada/SPM/Model2/ROI6/FC'
        path_to='H:/MRI/pnTTC/Prosociality_DC_Dr_Okada/SPM/Model2/ROI6/FC'
        siblings=[1,2,3,4]
        roi=6
        ###########

        for sibling in siblings:
            path_folder_to = path_to + '/Sibling' + str(sibling)
            if not os.path.exists(path_folder_to):
                os.makedirs(path_folder_to)
            with open(path_from+ '/id_sibling' + str(sibling) + '.txt', 'r') as list_id:
                list_id=list_id.readlines()
                list_id=[x.strip('\n') for x in list_id]
            for i in list_id:
                path_file_from=path_from + '/zROI' + str(roi) + 'FCMap_CSUB-' + str(i).zfill(5) + 'C-01.nii'
                if os.path.exists(path_file_from):
                    shutil.copy(path_file_from,path_folder_to)
                else:
                    print('File ' + path_file_from + ' does not exist.')           


#######################
# CONVERT PICKUP FILE #
#######################
# convert id according to list and pickup file

class ConvertPickupFile():
    def __init__(self):
        ############
        # Parameters
        path_from=''
        path_to=''
        file_id=''
        file_id_convert=''
        prefix=''
        suffix=''
        ############

        df_id = pd.read_csv(path_from + '/' + file_id)
        df_convert=pd.read_csv(path_from + '/' + file_id_convert)
        df_name=pd.merge(df_id,df_convert, how='left', on='ID_exam')
        for name in df_name['ID_MRI']:
            name_file=name + suffix
            file_from=path_from + '/' + name_file
            file_to=path_to + '/' + name_file
            shutil.copy(file_from, file_to)
            print('Copied file: ' + file_from)
        print('Finished file pick-up.')


############################
# COPY FOLDER NAME TO FILE #
############################

class Folder2File():
    def __init__(self):
        ############
        # Parameters
        path_from=''
        path_to=''
        ############
        list_file=glob.glob(path_from + '/*/output.nii')
        for path in list_file:
            name_subfolder=path.split('\\')[1]
            file_from=path_from + '/' + name_subfolder + '/output.nii'
            file_to=path_to + '/'+ name_subfolder + '.nii'
            shutil.copy(file_from, file_to)
            print('Converted file: '+ file_from + '.')
        print('Done.')
#############
# LIBRARIES #
#############

import os
import shutil
import glob
import numpy as numpy
import pandas as pd

# copy file name to folder name and put file into folder
# used for dparsf
#suffix = '.nii'
#path='J:/MRI/pnTTC/Prosociality_DC_Dr_Okada/Analysis/FunImg'
#path='J:/MRI/pnTTC/Prosociality_DC_Dr_Okada/Analysis/T1Img'

class FileIntoFolder():
    def __init__(self,path=path,suffix=suffix):
        self.path=path
        self.suffix=suffix
        self.list_dir = os.listdir(self.path)
        for f in self.list_dir:
            if f.endswith(self.suffix):
                name_folder = f[:-(len(self.suffix))]
                path_newfolder = self.path + '/' + name_folder + '/'
                path_oldfile = self.path + '/' + f
                if not os.path.exists(path_newfolder):
                    os.makedirs(path_newfolder)
                shutil.move(path_oldfile, path_newfolder)
                print('Moved ' + f + '.')


# used for sorting DPARSF-analyzed data into folder
path_from='H:/MRI/pnTTC/Prosociality_DC_Dr_Okada/SPM/FC_FunRawARCWSF'
path_to='H:/MRI/pnTTC/Prosociality_DC_Dr_Okada/SPM'
path_list_id='H:/MRI/pnTTC/Prosociality_DC_Dr_Okada/SPM/Script'

#list_rois=[1,2,3,4,5,6,7,8]
list_rois=[6]
#list_models=[1,2,3]
list_models=[2]

class PickupROIFile():
    def __init__(self,path_from=path_from,path_to=path_to,path_list_id=path_list_id,
                 list_rois=list_rois,list_models=list_models):
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


path_from='H:/MRI/pnTTC/Prosociality_DC_Dr_Okada/SPM/Model2/ROI6/FC'
path_to='H:/MRI/pnTTC/Prosociality_DC_Dr_Okada/SPM/Model2/ROI6/FC'
siblings=[1,2,3,4]
roi=6
class SortROIFile():
    def __init__(self,path_from=path_from):
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

# for general file selection

class Pickup_File():
    def __init__(self,path_from=path_from,path_to=path_to,path_id=path_id,file_id=file_id,file_roi=file_roi,prefix=prefix,suffix=suffix):
        file_id_open=open(path_id + '/' + file_id, 'r')
        file_id_open=file_id_open.readlines()
        self.id_list=[int(x.strip('\n')) for x in file_id_open]
        for id in self.id_list:
            name_file=prefix + str(id).zfill(5) + suffix
            file_from=path_from + '/' + name_file
            file_to=path_to + '/' + name_file
            shutil.copy(file_from, file_to)
            print('Copied file: ' + file_from)
        print('Finished file pick-up.')


class Convert_Pickup_File():
    def __init__(self,path_from=path_from,path_to=path_to,file_id=file_id,file_id_convert=file_id_convert,prefix=prefix,suffix=suffix):
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


class Convert_ID():
    def __init__(self,path_from=path_from,file_id=file_id,file_id_convert=file_id_convert):
        df_id = pd.read_csv(path_from + '/' + file_id)
        df_convert=pd.read_csv(path_from + '/' + file_id_convert)
        df_name=pd.merge(df_id,df_convert, how='left', on='ID_exam')
        df_name.to_csv(path_from + '/ID_converted.csv')


class Folder2File():
    def __init__(self,path_from=path_from,path_to=path_to):
        list_file=glob.glob(path_from + '/*/output.nii')
        for path in list_file:
            name_subfolder=path.split('\\')[1]
            file_from=path_from + '/' + name_subfolder + '/output.nii'
            file_to=path_to + '/'+ name_subfolder + '.nii'
            shutil.copy(file_from, file_to)
            print('Converted file: '+ file_from + '.')
        print('Done.')

print('End of file')

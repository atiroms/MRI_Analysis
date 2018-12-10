

#path = '/media/atiroms/MORITA_HDD4/MRI/COCORO/03_analysis2/FunImg'
#path = '/media/atiroms/MORITA_HDD4/MRI/COCORO/test2/FunImg'
#path='/media/atiroms/MORITA_HDD4/MRI/COCORO/test2/Results/FC_FunImgARWSCF'
#path = '/media/veracrypt1/MRI/COCORO/01_analysis/FunImg'

#path_from = '/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/13_nii.gz_unite'
#path_to = '/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/14_qc'


#prefix = 'ROI'
#prefix = 'CSUB-'

#suffix = '.nii.gz'
#suffix = 'C-02.nii.gz'

#file_id = 'pnTTC2_T1_QC.txt'

#path_from = '/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/13_nii.gz_unite'
#path_to = '/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/14_qc'

###
#path_from = 'D:/atiroms/MRI/pnTTC/pnTTC2_rsfMRI_C/CONN/15_conn/T1_14_qc'
#path_to = 'D:/atiroms/MRI/pnTTC/pnTTC2_rsfMRI_C/CONN/15_conn/T1'
#path_from = 'D:/atiroms/MRI/pnTTC/pnTTC2_rsfMRI_C/CONN/15_conn/rsfMRI_13_nii.gz_unite'
#path_to = 'D:/atiroms/MRI/pnTTC/pnTTC2_rsfMRI_C/CONN/15_conn/rsfMRI'
#prefix = 'CSUB-'
#suffix = 'C-02.nii.gz'
#file_id = 'ID_W2_T1QC_rsfMRIexist.txt'

###
#path_from='/media/veracrypt1/MRI/COCORO/02_analysis_f/Results/FC_FunImgARWSCF'
#path_from='/media/veracrypt1/MRI/COCORO/02_analysis_f/Results/FC_FunImgARWSglobalCF'
#path_to='/media/veracrypt1/MRI/COCORO/02_analysis_f/Results/zROI_FC_FunImgARWSCF'
#path_to='/media/veracrypt1/MRI/COCORO/02_analysis_f/Results/zROI_FC_FunImgARWSglobalCF'
#path_list_id='/media/veracrypt1/MRI/COCORO'

#rois = [37, 38, 41, 42, 71, 72, 73, 74, 75, 76, 77, 78]
#groups = [0,1,3]

###
#path_from = 'F:/MRI/Original_Data/Closed_Data/MR7/QC/tar.gz'
#path_to = 'F:/MRI/Original_Data/Closed_Data/MR7/QC/selected'
#prefix = ''
#suffix = '.tar.gz'
#file_id = 'ID.csv'
#file_id_convert = 'ID_convert.csv'

##
path_from = 'F:/MRI/Original_Data/Closed_Data/MR7/QC/nii'
path_to = 'F:/MRI/Original_Data/Closed_Data/MR7/QC/nii.gz'


import os
import shutil
import glob
import numpy as numpy
import pandas as pd

# copy file name to folder name and put file into folder
# used for dparsf
suffix = '.nii'
#path='J:/MRI/pnTTC/Prosociality_DC_Dr_Okada/Analysis/FunImg'
path='J:/MRI/pnTTC/Prosociality_DC_Dr_Okada/Analysis/T1Img'

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
'''
class PickupROIFile_Old():
    def __init__(self,path=path,prefix=prefix,rois=rois):
        self.path=path
        self.prefix=prefix
        for roi in rois:
            list_file=glob.glob(self.path + '/zROI' + str(roi) + 'FC*.nii')
            path_tofolder = self.path + '/zROI' + str(roi)
            if not os.path.exists(path_tofolder):
                os.makedirs(path_tofolder)
            for f in list_file:
                shutil.copy(f,path_tofolder)
            print('Picked ROI ' + str(roi) + ' files.')
'''

class PickupROIFile():
    def __init__(self,path_from=path_from,path_to=path_to,path_list_id=path_list_id,rois=rois,groups=groups):
        for roi in rois:
            for group in groups:
                path_tofolder = path_to + '/zROI' + str(roi) + '_' + str(group)
                if not os.path.exists(path_tofolder):
                    os.makedirs(path_tofolder)
                file=open(path_list_id+ '/ID_group' + str(group) + '.txt', 'r')
                file=file.readlines()
                list_id=[x.strip('\n') for x in file]
                for i in list_id:
                    path_file_from=path_from + '/zROI' + str(roi) + 'FCMap_F_' + i + '.nii'
                    if os.path.exists(path_file_from):
                        shutil.copy(path_file_from,path_tofolder)
                print('Copied ROI ' + str(roi) + ', Group ' + str(group) + ' files.')
        print('Done.')


# for general use
#path_from='G:/MRI/pnTTC/pnTTC1_rsfMRI_C/CONN/04_nii'
path_from='G:/MRI/pnTTC/pnTTC1_T1_C/FS/04_nii'
#path_to='G:/MRI/pnTTC/pnTTC1_rsfMRI_C/CONN/FunRaw'
path_to='G:/MRI/pnTTC/pnTTC1_T1_C/FS/T1Img'
file_id='birthorder_QC.txt'
prefix = 'CSUB-'
suffix = 'C-01.nii'

class Pickup_File():
    def __init__(self,path_from=path_from,path_to=path_to,file_id=file_id,prefix=prefix,suffix=suffix):
        file=open(path_from + '/' + file_id, 'r')
        file=file.readlines()
        self.id_list=[int(x.strip('\n')) for x in file]
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

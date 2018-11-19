

#path = '/media/atiroms/MORITA_HDD4/MRI/COCORO/03_analysis2/FunImg'
#path = '/media/atiroms/MORITA_HDD4/MRI/COCORO/test2/FunImg'
#path='/media/atiroms/MORITA_HDD4/MRI/COCORO/test2/Results/FC_FunImgARWSCF'
#path = '/media/veracrypt1/MRI/COCORO/01_analysis/FunImg'

#path_from = '/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/13_nii.gz_unite'
#path_to = '/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/14_qc'

#rois = [37, 38, 41, 42, 71, 72, 73, 74, 75, 76, 77, 78]
#prefix = 'ROI'
#prefix = 'CSUB-'

#suffix = '.nii.gz'
#suffix = 'C-02.nii.gz'

#file_id = 'pnTTC2_T1_QC.txt'

#path_from = '/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/13_nii.gz_unite'
#path_to = '/media/veracrypt1/MRI/pnTTC/pnTTC2_T1_C/FS/14_qc'

#path_from = 'D:/atiroms/MRI/pnTTC/pnTTC2_rsfMRI_C/CONN/15_conn/T1_14_qc'
#path_to = 'D:/atiroms/MRI/pnTTC/pnTTC2_rsfMRI_C/CONN/15_conn/T1'
path_from = 'D:/atiroms/MRI/pnTTC/pnTTC2_rsfMRI_C/CONN/15_conn/rsfMRI_13_nii.gz_unite'
path_to = 'D:/atiroms/MRI/pnTTC/pnTTC2_rsfMRI_C/CONN/15_conn/rsfMRI'
prefix = 'CSUB-'
suffix = 'C-02.nii.gz'
file_id = 'ID_W2_T1QC_rsfMRIexist.txt'

import os
import shutil
import glob

# copy file name to folder name and put file into folder
# used for dparsf
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
class PickupROIFile():
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



# for general use
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


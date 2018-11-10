

path = '/media/atiroms/MORITA_HDD4/MRI/COCORO/03_analysis2/FunImg'
#path = '/media/atiroms/MORITA_HDD4/MRI/COCORO/test2/FunImg'
#path='/media/atiroms/MORITA_HDD4/MRI/COCORO/test2/Results/FC_FunImgARWSCF'

rois = [37, 38, 41, 42, 71, 72, 73, 74, 75, 76, 77, 78]
prefix = 'ROI'

suffix = '.nii.gz'

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


class PickupFile():
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


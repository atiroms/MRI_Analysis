

path = '/media/atiroms/MORITA_HDD4/MRI/COCORO/Analysis/FunImg'
suffix = '.nii.gz'

import os
import shutil

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
            
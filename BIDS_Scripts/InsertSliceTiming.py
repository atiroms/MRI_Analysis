
##############
# Parameters #
##############

TR=2.5
n_slices=40
#path_exp='/media/veracrypt1/MRI/pnTTC/BIDS/Nifti'
path_exp='C:/Users/atiro/Dropbox/BIDS/test_3sub/Nifti2_SliceTiming'
sessions=['ses-01','ses-02']


#############
# LIBRARIES #
#############

import os
import json


#############
# MAIN CODE #
#############

class InsertSliceTiming():
    def __init__(self,TR=TR,n_slices=n_slices,path_exp=path_exp,sessions=sessions):
        self.path_exp=path_exp
        self.sessions=sessions
        list_dir_all = os.listdir(self.path_exp)
        self.list_dir_sub=[]
        self.list_dir_func=[]
        self.list_slicetiming=[]
        for i in range(n_slices):
            self.list_slicetiming.append(i*TR/n_slices)
        for dir_sub in list_dir_all:
            if dir_sub.startswith('sub-'):
                self.list_dir_sub.append(dir_sub)
                for session in self.sessions:
                    dir_func=self.path_exp +'/' + dir_sub + '/' + session + '/func'
                    if os.path.exists(dir_func):
                        self.list_dir_func.append(dir_func)
                        filename_json = dir_sub + '_' + session + '_task-rest_bold.json'
                        with open(dir_func + '/' + filename_json) as file_json_input:  
                            data = json.load(file_json_input)
                        data['SliceTiming']=self.list_slicetiming
                        with open(dir_func + '/' + filename_json, 'w') as file_json_output:  
                            json.dump(data, file_json_output,indent=2, sort_keys=True)
                        print('Added SliceTiming data to ' + filename_json + '.')

##############
# PARAMETERS #
##############

path_exp =''
file_mask = ''
dir_input = ''
dir_output = ''

#############
# LIBRARIES #
#############

import numpy as np
import os
from nilearn.masking import apply_mask


#############
# MAIN CODE #
#############

class ExtractAllVoxel():
    def __init__(self,path_exp=path_exp,file_mask=file_mask, dir_input=dir_input, dir_output=dir_output):
        list_input = os.listdir(path_exp + '/' + dir_input)
        path_mask = path_exp + '/' + file_mask
        for file_input in list_input:
            path_input=path_exp + '/' + dir_input + '/' + file_input
            data_masked = apply_mask(path_input, path_mask)
            file_output=path_input.replace('.nii.gz','.csv')
            path_output=path_exp + '/' + dir_output + '/' + file_output
            np.savetxt(path_output, data_masked, delimiter=",")
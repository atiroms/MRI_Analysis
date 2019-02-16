##############
# PARAMETERS #
##############

path_exp='D:/atiroms/MRI/pnTTC/BIDS/practice'
file_img='sub-00014_ses-01_task-rest_space-MNI152NLin2009cAsym_boldref.nii.gz'

#############
# LIBRARIES #
#############

from nilearn import plotting
import os


#############

path_img=os.path.join(path_exp,file_img)
plotting.plot_glass_brain(path_img)
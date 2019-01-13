
path_from='H:/MRI/pnTTC/Prosociality_DC_Dr_Okada/test'
file_from='0vs1.mat'

#############
# LIBRARIES #
#############

import scipy.io
import os

#############
# MAIN CODE #
#############

path_file_input=os.path.join(path_from, file_from)
spm_input = scipy.io.loadmat(path_file_input)

x=spm_input['matlabbatch'][0][0][0][0][0][0][0][0][0][0][0][0][0]

type(x)

x[0]

scipy.io.savemat('newmat.mat',spm_input)

x=spm_input['matlabbatch'][0,0]['spm'][0,0]['stats'][0,0]['factorial_design']



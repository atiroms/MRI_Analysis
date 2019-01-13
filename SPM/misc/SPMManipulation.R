
#### Parameters ####

path.input='H:/MRI/pnTTC/Prosociality_DC_Dr_Okada/test'
file.input='0vs1.mat'


#### Libraries ####

library(R.matlab)


#### Main Code ####

data.input<-readMat(file.path(path.input,file.input))

writeMat(file.path(path.input,'new2.mat'),data.input['matlabbatch'])

a<-data.input['matlabbatch']


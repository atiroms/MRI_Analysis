library(voxel)
library(mgcv)
library(oro.nifti)


#path_exp<-"C:/Users/atiro/Dropbox/Temp/5sub_corr"
path_exp<-"D:/atiroms/Dropbox/Temp/5sub_corr"
#list_path_src<-NULL
list_path_src<-list.files(file.path(path_exp,'input'))
list_path_src<-file.path(path_exp,'input',list_path_src)

#img_mask<-file.path(path_exp,'mask','mni_icbm152_brainmask_tal_nlin_asym_09c.nii')
img_mask<-file.path(path_exp,'mask','mni_icbm152_brain_resample_mask_tal_nlin_asym_09c.nii')
img_mask<-readNIfTI(img_mask)
img_mask<-img_mask>0
#img_mask<-file.path(path_exp,'mask','output.nii')
df_subjdata<-read.csv(file.path(path_exp,'df','data.csv'))

dst<-gamNIfTI(list_path_src,
              mask=img_mask,
              formula='~ value',
              subjData=df_subjdata)



#vgamParam(list_path_src,
#          mask=img_mask,
#          formula='~value',
#          subjData=df_subjdata)

#image <- mergeNiftis(inputPaths = list_path_src, direction = "t", outfile = NULL)
#mask <- oro.nifti::readNIfTI(fname=img_mask)
#imageMat <- ts2matrix(image, mask)
#label <- sort(as.numeric(unique(matrix(img_mask@.Data))))

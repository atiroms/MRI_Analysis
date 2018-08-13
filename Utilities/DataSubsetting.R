
# R script to subset data


#### Parameters ####

#working_dir<-"D:/MRI/Statistics/Functional_CONN_HO"
#working_dir<-"D:/MRI/Statistics/ConnectionClinical"
#working_dir<-"D:/MRI/Statistics/Structural_FS"
#working_dir<-"G:/MRI/Statistics/Functional_CONN_Power"
working_dir<-"G:/MRI/Statistics/Connection"

#input_file<-"FC_r.csv"
#input_file<-"Freesurfer_Volume.csv"
#input_file<-"Freesurfer_CorticalThickness.csv"
#input_file<-"Freesurfer_CorticalArea.csv"
#input_file<-"CONN_DK_FC_Subset182.csv"
#input_file<-"182_FS_DK_Area_JKPV.csv"
#input_file<-"182_FS_DK_Area_JKZ.csv"
input_file<-"DK_Area_JKPE.csv"

#subset_file<-"T1_rsfMRI_TS_male_list.csv"
subset_file<-"T1_rsfMRI_TS_female_list.csv"
#subset_file<-"T1_TS_male_list.csv"
#subset_file<-"T1_TS_female_list.csv"

#output_file<-"DK_FC_r_TSexist_male.csv"
#output_file<-"DK_FC_r_TSexist_female.csv"
#output_file<-"FS_Volume_TSexist_male.csv"
#output_file<-"FS_Volume_TSexist_female.csv"
#output_file<-"FS_Thickness_TSexist_male.csv"
#output_file<-"FS_Thickness_TSexist_female.csv"
#output_file<-"FS_Area_TSexist_male.csv"
#output_file<-"FS_Area_TSexist_female.csv"
#output_file<-"Power_FC_r_TSexist_male.csv"
#output_file<-"Power_FC_r_TSexist_female.csv"
#output_file<-"DK_Area_JKPV_TSexist_male.csv"
#output_file<-"DK_Area_JKPV_TSexist_female.csv"
#output_file<-"DK_Area_JKZ_TSexist_male.csv"
#output_file<-"DK_Area_JKZ_TSexist_female.csv"
#output_file<-"DK_Area_JKPE_TSexist_male.csv"
output_file<-"DK_Area_JKPE_TSexist_female.csv"

#### Data loading ####
input_data <-read.csv(file.path(working_dir, input_file), check.names=F)
#subset_data<-read.csv(file.path(working_dir, subset_file), check.names=F)
subset_data<-read.csv(file.path(working_dir, subset_file))
subset_data<-subset_data$ID_pnTTC
output_data<-data.frame(matrix(nrow=0, ncol=ncol(input_data)))
for (i in subset_data){
  output_data<-rbind(output_data,input_data[which(input_data$ID_pnTTC==i),])
}


write.csv(output_data, file.path(working_dir,output_file),row.names = F)
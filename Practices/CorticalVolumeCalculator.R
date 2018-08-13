##

working_dir <- "G:/MRI/Statistics/CorticalStructure"
thickness_file <- "Freesurfer_CorticalThickness02.csv"
area_file <- "Freesurfer_CorticalArea01.csv"
volume_file<-"Freesurfer_CorticalVolume01.csv"


setwd(working_dir)
thickness_data <-read.csv(thickness_file)
area_data <-read.csv(area_file)

volume_data<-thickness_data

for (i in 2:ncol(thickness_data)){
  for (j in 1:nrow(thickness_data)){
    volume_data[j,i]<-thickness_data[j,i]*area_data[j,i]
  }
}

write.csv(volume_data, file.path(working_dir,volume_file))
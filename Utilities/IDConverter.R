#

id_dir <- "G:/MRI/Info"
id_filename<-"ID6.csv"
inputidlabel<-"ID_T1QC_rsfMRIexist"
outputidlabel<-"ID_pnTTC"

id_list<-read.csv(file.path(id_dir,id_filename))
idlist_colname<-data.frame(label=colnames(id_list))

ConvertID<-function(input,inputidlabel,outputidlabel){
  input_col<-id_list[,which(idlist_colname$label==inputidlabel)]
  output_col<-id_list[,which(idlist_colname$label==outputidlabel)]
  output<-data.frame(matrix(nrow=length(input),ncol=1))
  for (i in 1:length(input)){
    input_value<-input[[i]]
    output_row<-which(input_col==input_value)
    output[i,1]<-output_col[output_row]
  }
  output<-as.numeric(output[,1])
  return(output)
}
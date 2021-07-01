library(easypackages)
libraries('data.table',"dplyr")

path_src<-'D:/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/common/clinical_data_manipulation/source'
path_parent<-'D:/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/common/clinical_data_manipulation'
file_ib<-'191027TTC_itembank_labelling.csv'
file_spss<-'spss_data.csv'

df_ib<-data.frame(fread(file.path(path_src,file_ib),encoding='UTF-8'))
df_spss<-data.frame(fread(file.path(path_src,file_spss)))


df_item<-data.frame(item_id=colnames(df_spss))

df_item<-dplyr::left_join(df_item,df_ib,by=c("item_id"=colnames(df_ib)[2]))

fwrite(df_item,file.path(path_src,'items01.csv'))

# manually delete duplicate rows, change column names

df_item<-data.frame(fread(file.path(path_parent,'items02.csv'),encoding='UTF-8'))

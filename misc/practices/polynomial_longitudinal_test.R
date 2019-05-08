
path_exp <- "Dropbox/MRI/pnTTC/Puberty/Stats/T1w_FS"
dir_in <-"01_extract"
dir_out <-"04_gamm_test"
file_input<-"fs_measure.csv"

list_wave <- c(1,2)

#list_measure <-c("volume","thickness","area")
list_measure <-"volume"

#list_str_group<-c("cortex","subcortex","white matter","global","misc")
list_str_group<-"subcortex"
#list_str_group<-c("global","misc")
#list_str_group<-c("cortex","subcortex","global")

#list_covar<-list("tanner"=list("1"="W1_Tanner_Max",
#                               "2"="W2_Tanner_Max",
#                               "label"="Tanner stage"),
#                 "age"=list("1"="W1_Age_at_MRI",
#                            "2"="W2_Age_at_MRI",
#                            "label"="Age"))
#list_covar<-list("tanner"=list("1"="W1_Tanner_Full",
#                               "2"="W2_Tanner_Full",
#                               "label"="Tanner stage"),
#                 "age"=list("1"="W1_Age_at_MRI",
#                            "2"="W2_Age_at_MRI",
#                            "label"="Age"))
list_covar<-list("tanner"=list("1"="W1_Tanner_Max",
                               "2"="W2_Tanner_Max",
                               "label"="Tanner stage"),
                 "age"=list("1"="W1_Age_at_MRI",
                            "2"="W2_Age_at_MRI",
                            "label"="Age"),
                 "sex"=list("1"="Sex",
                            "2"="Sex",
                            "label"="Sex"))


#subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
#                             list("key"="W1_T1QC_new_mild","value"=1),
#                             list("key"="Sex","value"=1)),
#                    "2"=list(list("key"="W2_T1QC","value"=1),
#                             list("key"="W2_T1QC_new_mild","value"=1),
#                             list("key"="Sex","value"=1)))

subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
                             list("key"="W1_T1QC_new_mild","value"=1)),
                    "2"=list(list("key"="W2_T1QC","value"=1),
                             list("key"="W2_T1QC_new_mild","value"=1)))

#subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
#                             list("key"="Sex","value"=2)),
#                    "2"=list(list("key"="W2_T1QC","value"=1),
#                             list("key"="Sex","value"=2)))

#subset_subj <- list("1"=list(list("key"="W1_T1QC_new_mild","value"=1),
#                             list("key"="Sex","value"=2)),
#                    "2"=list(list("key"="W2_T1QC_new_mild","value"=1),
#                             list("key"="Sex","value"=2)))

#str_mod <- "value ~ s(age,k=3) + s(tanner,k=3) + s(ID_pnTTC,bs='re')"
#str_mod <- "value ~ s(age) + s(tanner,k=3) + s(ID_pnTTC,bs='re')"
#str_mod <- "value ~ age + tanner + s(ID_pnTTC,bs='re')"
#str_mod <- "value ~ age*tanner + s(ID_pnTTC,bs='re')"
#str_mod <- "value ~ age + tanner:sex + s(ID_pnTTC,bs='re')"
#list_mod <- list("a+ts"=
#                   "value ~ age + tanner:sex + s(ID_pnTTC,bs='re')",
#                 "a+ts+ats"=
#                   "value ~ age + tanner:sex + age:tanner:sex + s(ID_pnTTC,bs='re')")
#list_mod <- list("s+a+ts"=
#                   "value ~ sex + age + tanner:sex + s(ID_pnTTC,bs='re')",
#                 "s+a+ts+ats"=
#                   "value ~ sex + age + tanner:sex + age:tanner:sex + s(ID_pnTTC,bs='re')")
list_mod <- list("s+a+ts"=
                   "value ~ sex + age + poly(tanner,2):sex + s(ID_pnTTC,bs='re')",
                 "s+a+ts+ats"=
                   "value ~ sex + age + poly(tanner,2):sex + age:poly(tanner,2):sex + s(ID_pnTTC,bs='re')")

list_graph <- list("s+a+ts"=
                     list("a"=list("title"="Age effect",
                                   "x_axis"="age",
                                   "smooth"=list("Male"=list("fix"=list("sex"=1),
                                                             "color"="steelblue2","alpha"=1,"ribbon"=T),
                                                 "Female"=list("fix"=list("sex"=2),
                                                               "color"="lightcoral","alpha"=1,"ribbon"=T)),
                                   "point"=list("Male"=list("subset"=list("sex"=1),
                                                            "color"="steelblue2","alpha"=1),
                                                "Female"=list("subset"=list("sex"=2),
                                                              "color"="lightcoral","alpha"=1))),
                          "ts"=list("title"="Tanner stage effect",
                                    "x_axis"="tanner",
                                    "smooth"=list("Male"=list("fix"=list("sex"=1),
                                                              "color"="steelblue2","alpha"=1,"ribbon"=T),
                                                  "Female"=list("fix"=list("sex"=2),
                                                                "color"="lightcoral","alpha"=1,"ribbon"=T)),
                                    "point"=list("Male"=list("subset"=list("sex"=1),
                                                             "color"="steelblue2","alpha"=1),
                                                 "Female"=list("subset"=list("sex"=2),
                                                               "color"="lightcoral","alpha"=1)))),
                   "s+a+ts+ats"=
                     list("a"=list("title"="Age effect",
                                   "x_axis"="age",
                                   "smooth"=list("Male"=list("fix"=list("sex"=1),
                                                             "color"="steelblue2","alpha"=1,"ribbon"=T),
                                                 "Female"=list("fix"=list("sex"=2),
                                                               "color"="lightcoral","alpha"=1,"ribbon"=T)),
                                   "point"=list("Male"=list("subset"=list("sex"=1),
                                                            "color"="steelblue2","alpha"=1),
                                                "Female"=list("subset"=list("sex"=2),
                                                              "color"="lightcoral","alpha"=1))),
                          "ts"=list("title"="Tanner stage effect",
                                    "x_axis"="tanner",
                                    "smooth"=list("Male"=list("fix"=list("sex"=1),
                                                              "color"="steelblue2","alpha"=1,"ribbon"=T),
                                                  "Female"=list("fix"=list("sex"=2),
                                                                "color"="lightcoral","alpha"=1,"ribbon"=T)),
                                    "point"=list("Male"=list("subset"=list("sex"=1),
                                                             "color"="steelblue2","alpha"=1),
                                                 "Female"=list("subset"=list("sex"=2),
                                                               "color"="lightcoral","alpha"=1))),
                          "ats"=list("title"="Age-Tanner stage interaction",
                                     "x_axis"="age",
                                     "smooth"=list("Male TS = 1"=list("fix"=list("sex"=1,"tanner"=1),
                                                                      "color"="Steelblue2","alpha"=0.4,"ribbon"=F),
                                                   "Male TS = 3"=list("fix"=list("sex"=1,"tanner"=3),
                                                                      "color"="steelblue2","alpha"=0.7,"ribbon"=F),
                                                   "Male TS = 5"=list("fix"=list("sex"=1,"tanner"=5),
                                                                      "color"="steelblue2","alpha"=1,"ribbon"=F),
                                                   "Female TS = 1"=list("fix"=list("sex"=2,"tanner"=1),
                                                                        "color"="lightcoral","alpha"=0.4,"ribbon"=F),
                                                   "Female TS = 3"=list("fix"=list("sex"=2,"tanner"=3),
                                                                        "color"="lightcoral","alpha"=0.7,"ribbon"=F),
                                                   "Female TS = 5"=list("fix"=list("sex"=2,"tanner"=5),
                                                                        "color"="lightcoral","alpha"=1,"ribbon"=F)),
                                     "point"=list("Male"=list("subset"=list("sex"=1),
                                                              "color"="steelblue2","alpha"=1),
                                                  "Female"=list("subset"=list("sex"=2),
                                                                "color"="lightcoral","alpha"=1)))))


#**************************************************
# Libraries =======================================
#**************************************************
library(Hmisc)
library(tidyverse)
library(tidyr)
library(mgcv)
library(dplyr)
library(ggplot2)


#**************************************************
# Create path list ================================
#**************************************************
func_path<-function(list_path_root = c("D:/atiroms","C:/Users/atiro","/home/atiroms","C:/Users/NICT_WS"),
                    path_exp_=path_exp,
                    dir_in_=dir_in,
                    dir_out_=dir_out){
  path_root<-NA
  for(p in list_path_root){
    if(file.exists(p)){
      path_root<-p
    }
  }
  if(is.na(path_root)){
    print("Error: root path could not be found.")
  }
  path_script <- file.path(path_root,"GitHub/MRI_Analysis")
  path_common <- file.path(path_root,"DropBox/MRI/pnTTC/Puberty/Stats/CommonData")
  path_in     <- file.path(path_root,path_exp_,dir_in_)
  path_out    <- file.path(path_root,path_exp_,dir_out_)
  output <- list("script"=path_script,"input"=path_in,"output"=path_out,
                 "common"=path_common,"dir_in"=dir_in_,"dir_out"=dir_out_)
  return(output)
}

paths<-func_path()


#**************************************************
# Original library ================================
#**************************************************
source(file.path(paths$script,"functionality/function.R"))
#source(file.path(paths$script,"functionality/glm_function.R"))
source(file.path(paths$script,"functionality/graph.R"))
#source(file.path(script_dir,"Functionalities/LI_Functions.R"))


#**************************************************
# GAMM of structural measures =====================
#**************************************************

paths_=paths
subset_subj_=subset_subj
list_covar_=list_covar
file_input_=file_input
list_wave_=list_wave
list_measure_=list_measure
list_str_group_=list_str_group
list_mod_=list_mod
list_graph_=list_graph


print("Starting gamm_str().")
nullobj<-func_createdirs(paths_,copy_log=T)

# Load and subset clinical data according to specified subsetting condition and covariate availability
print('Loading clinical data.')
df_clin<-func_clinical_data_long(paths_,list_wave_)
data_subset_clin<-func_subset_clin(df_clin,
                                   list_wave_,list_measure_,subset_subj_,
                                   list_covar_,
                                   rem_na_clin=T)
df_clin_subset<-data_subset_clin$df_clin

# Load structural data
print('Loading structural data.')
df_str<-read.csv(file.path(paths_$input,"output",file_input_))
df_str$value[which(is.nan(df_str$value))]<-0
dict_roi<-func_dict_roi(paths_)
df_str_subset<-func_subset_str(df_str,
                               list_measure_,
                               list_str_group_,
                               dict_roi)

# Join clinical and structural data frames
print('Joining clinical and structural data.')
df_join<-inner_join(df_str_subset,df_clin_subset,by=c('ID_pnTTC','wave'))
for (key in c('ID_pnTTC','wave','sex','measure')){
  if (key %in% colnames(df_join)){
    df_join[,key]<-as.factor(df_join[,key])
  }
}

# Calculate GAMM
print('Calculating GAMM.')
df_out_term<-data.frame(matrix(nrow=0,ncol=9))
colnames(df_out_term)<-c("measure","roi","label_roi","group_roi","model","term","F","t","p")
df_out_model<-data.frame(matrix(nrow=0,ncol=6))
colnames(df_out_model)<-c("measure","roi","label_roi","group_roi","F","p")

measure<-list_measure_[1]

print(paste('Calculating measurements of ',measure,sep=''))
df_join_measure<-df_join[df_join$measure==measure,]
list_roi<-as.character(unique(df_join_measure$roi))
list_roi<-list_roi[order(list_roi)]

roi<-list_roi[1]

label_roi<-as.character(dict_roi[dict_roi$id==roi,'label'])
group_roi<-as.character(dict_roi[dict_roi$id==roi,'group'])
print(paste('Calculating ',roi,' = ',label_roi,sep=''))
df_join_measure_roi<-df_join_measure[df_join_measure$roi==roi,]
mod_gamm<-list()

mod<-names(list_mod_)[1]

mod_gamm[[mod]]<-gam(as.formula(list_mod_[[mod]]),data=df_join_measure_roi)
p_table<-summary.gam(mod_gamm[[mod]])$p.table
s_table<-summary.gam(mod_gamm[[mod]])$s.table
df_out_term_add<-rbind(data.frame(measure=measure,roi=roi,label_roi=label_roi,group_roi=group_roi,model=mod,
                                  term=rownames(p_table),F=NA,t=p_table[,'t value'],p=p_table[,'Pr(>|t|)']),
                       data.frame(measure=measure,roi=roi,label_roi=label_roi,group_roi=group_roi,model=mod,
                                  term=rownames(s_table),F=s_table[,'F'],t=NA,p=s_table[,'p-value']))
df_out_term<-rbind(df_out_term,df_out_term_add)
list_graph_mod<-list_graph_[[mod]]

idx_graph<-names(list_graph_mod)[2]

mod_gamm<-mod_gamm[[mod]]
spec_graph<-list_graph_mod[[idx_graph]]


#df_src <- mod_gamm$model

df_src <- df_join_measure_roi[names(mod_gamm$var.summary)]

plot<-ggplot()

name_smooth<-names(spec_graph[["smooth"]])[1]

spec_smooth<-spec_graph[["smooth"]][[name_smooth]]
df_smooth <- data.frame(x = seq(min(df_src[spec_graph[["x_axis"]]]),
                                max(df_src[spec_graph[["x_axis"]]]),
                                length.out=200))
names(df_smooth) <- spec_graph[["x_axis"]]
for (i in names(df_src)) {
  if (i != spec_graph[["x_axis"]]) {
    if (any(class(df_src[i][,1])[1] == c("numeric", "integer","boolean"))) {
      df_smooth[, dim(df_smooth)[2] + 1] <- mean(df_src[i][,1])
      names(df_smooth)[dim(df_smooth)[2]] <- i
    }
    else if (any(class(df_src[i][,1])[1] == c("character", "factor","ordered"))) {
      df_smooth[, dim(df_smooth)[2] + 1] <- df_src[i][1,1]
      names(df_smooth)[dim(df_smooth)[2]] <- i
    }
  }
}
if (!is.null(spec_smooth[["fix"]])){
  for (var in names(spec_smooth[["fix"]])){
    df_smooth[[var]]<-spec_smooth[["fix"]][[var]]
  }
}
#df_smooth <- cbind(df_smooth,
#                   as.data.frame(predict.gam(mod_gamm, df_smooth, se.fit = TRUE)))
df_smooth <- cbind(df_smooth,
                   as.data.frame(predict.gam(mod_gamm, df_smooth, exclude="s(ID_pnTTC)",se.fit = TRUE)))
plot <- (plot
         + geom_line(aes(x=!!df_smooth[,1],y=!!df_smooth[["fit"]]),
                     color=spec_smooth[["color"]],size=0.5,alpha=spec_smooth[["alpha"]]))


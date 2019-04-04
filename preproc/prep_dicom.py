##################################################
# Libraries
##################################################

import os
import shutil
import pandas as pd
import csv
#import nilearn.image as nl_image
#import json
import numpy as np
import pydicom
#import datetime
#import gzip
import tarfile
import glob


#def _copyfileobj_patched(fsrc, fdst, length=16*1024*1024):
def _copyfileobj_patched(fsrc, fdst, length=1024*1024*1024):
    """Patches shutil method to hugely improve copy speed"""
    while 1:
        buf = fsrc.read(length)
        if not buf:
            break
        fdst.write(buf)
shutil.copyfileobj = _copyfileobj_patched


##################################################
# make tar.gz file of sequence subfolders
##################################################
# used to manipulate pnTTC raw data

class TarGz():
    def __init__(self,
        #path_in='/media/veracrypt2/MRI/pnTTC/Raw/HUMAN-01-ANON_test',
        #path_out='/media/veracrypt1/MRI/pnTTC/Raw/HUMAN-01-ANON_test',
        #path_in='/media/veracrypt2/MRI/pnTTC/Raw/HUMAN-01-ANON',
        #path_out='/media/veracrypt1/MRI/pnTTC/Raw/HUMAN-01-ANON_zip',
        path_in='/Volumes/MRI_Ext1/smorita/MRI/pnTTC/Raw/HUMAN-01-ANON',
        path_out='/Volumes/MRI_Ext1/smorita/MRI/pnTTC/Raw/HUMAN-01-ANON_zip',
        type_subj=['C-01','C-02','M-01','P-01']
        ):

        print('Starting to pickup subjects and compressing files.')
        list_dir = os.listdir(path_in)
        list_dir_slctd=[]
        for t_s in type_subj:
            list_dir_slctd_add=[f for f in list_dir if t_s in f]
            list_dir_slctd_add.sort()
            list_dir_slctd =list_dir_slctd+list_dir_slctd_add
        print('Number of subjects / studies: ' + str(len(list_dir_slctd)))

        for dir_subj in list_dir_slctd:
            path_dir_subj_out=os.path.join(path_out,dir_subj)
            os.mkdir(path_dir_subj_out)
            list_dir_study=os.listdir(os.path.join(path_in,dir_subj))
            list_dir_study=[d for d in list_dir_study if d!='.DS_Store']
            list_dir_seq=os.listdir(os.path.join(path_in,dir_subj,list_dir_study[0]))
            for dir_seq in list_dir_seq:
                path_dir_in=os.path.join(path_in, dir_subj, list_dir_study[0], dir_seq)
                path_file_out=os.path.join(path_dir_subj_out,dir_seq+'.tar.gz')
                with tarfile.open(path_file_out, "w:gz") as tar:
                    tar.add(path_dir_in, arcname=os.path.basename(path_dir_in))
            print('Finished compressiong ' + dir_subj)
        
        print('Finished pickup and compressing files.')


##################################################
# Extract tar.gz file of sequence subfolders
##################################################
# used as HeuDiConv preparation

class UntarGz():
    def __init__(self,
        path_in='/media/veracrypt2/MRI/pnTTC/Raw/HUMAN-01-ANON_zip',
        path_out='/media/veracrypt1/MRI/pnTTC/Preproc/35_c2_dcm_slctd',
        #list_type_subj=['C-01','C-02','M-01'],
        list_type_subj=['C-02'],
        list_type_sequence=['+MPRAGE_CBSN -','+rsfMRI_SBPRS -','+rsfMRI -','+Fieldmap_SBPRS -','+Fieldmap -']
        ):

        print('Starting to Untar files')

        # Create experiment folder
        print('Starting to create experiment folder.')
        list_paths_mkdir=[]
        list_paths_mkdir.append(path_out)
        list_paths_mkdir.append(os.path.join(path_out,'output'))
        for type_subj in list_type_subj:
            list_paths_mkdir.append(os.path.join(path_out,'output',type_subj))
        for p in list_paths_mkdir:
            if not os.path.exists(p):
                os.makedirs(p)
        print('Finished creating experiment folder.')
        
        for type_subj in list_type_subj:
            print('Extracting '+type_subj)
            list_dir_subj = os.listdir(os.path.join(path_in,type_subj))
            list_dir_subj.sort()
            print('Number of subject folders: ' + str(len(list_dir_subj)))

            for dir_subj in list_dir_subj:
                id_subj=int(dir_subj[5:10])
                dir_subj_out='CSUB-'+str(id_subj).zfill(5)+type_subj
                path_dir_subj_out=os.path.join(path_out,'output',type_subj,dir_subj_out)
                if not os.path.exists(path_dir_subj_out):
                    os.makedirs(path_dir_subj_out)
                for type_sequence in list_type_sequence:
                    list_dir_sequence=glob.glob(os.path.join(path_in,type_subj,dir_subj)+'/'+type_sequence+'*')
                    if len(list_dir_sequence)>0:
                        if len(list_dir_sequence)>1:
                            print('Multiple sequence folders for subject: '+dir_subj+', sequence: '+type_sequence)
                        for i in range(len(list_dir_sequence)):
                            file_tar=tarfile.open(list_dir_sequence[i])
                            file_tar.extractall(path=path_dir_subj_out)
                            file_tar.close()
                    #elif len(list_dir_sequence)==0:
                    #    print('No sequence folders for subject: '+dir_subj+', sequence: '+type_sequence)
        
        print('Finished Untarring files')


##################################################
# Extract DICOM header metadata
##################################################
# used to extract scan date information

class ExtractDcmHeader():
    def __init__(self,
        #path_file_in='P:/MRI/pnTTC/Preproc/00_dicom_ses12_exist/pnTTC2_T1/CSUB-00003C-02/IM-0001-0001-0001.dcm'
        #path_file_in='/media/veracrypt1/MRI/pnTTC/Preproc/test_1sub/30_heudiconv/input/CSUB-00014C-01/CSUB-00014C-01/Csub/+Fieldmap_SBPRS - 301/IM-0001-0001-0001.dcm',
        #path_file_in='/media/veracrypt1/MRI/pnTTC/Preproc/test_1sub/30_heudiconv/input/CSUB-00014C-01/CSUB-00014C-01/Csub/+Fieldmap_SBPRS - 301/IM-0001-0002-0001.dcm',
        #path_file_out='/media/veracrypt1/MRI/pnTTC/Preproc/test_1sub/31_dcm2nii/output/header/output.csv',
        path_file_in='/media/veracrypt1/MRI/pnTTC/Preproc/test_1sub/30_heudiconv/input/CSUB-00014C-01/CSUB-00014C-01/Csub/+rsfMRI_SBPRS - 401/IM-0001-0001-0001.dcm',
        path_file_out='/media/veracrypt1/MRI/pnTTC/Preproc/test_1sub/31_dcm2nii/output/rsfMRI_header/header.csv',
        output_file=True
        ):

        img_in=pydicom.read_file(path_file_in)
        df_header=pd.DataFrame(columns=['tag','name','value'])
        for element in img_in:
            df_header=df_header.append(pd.Series([str(element.tag),element.name,element.repval],index=df_header.columns),ignore_index=True)
        if output_file:
            df_header.to_csv(path_file_out)
        self.output=df_header
        #print('DICOM element count: '+ str(len(df_header)))
        #print('All done.')


class ExtractAll():
    def __init__(self,
        path_in='/media/veracrypt1/MRI/pnTTC/Preproc/test_1sub/30_heudiconv/input/CSUB-00014C-01/CSUB-00014C-01/Csub/+Fieldmap_SBPRS - 301',
        path_file_out='/media/veracrypt1/MRI/pnTTC/Preproc/test_1sub/31_dcm2nii/output/header/header.csv',
        n_keys=353
        ):
        list_file = os.listdir(path_in)
        list_file.sort()
        #label_col=['file']+[str(num) for num in np.arange(n_keys)]
        #df_header_mlt=pd.DataFrame(columns=label_col)
        df_header_mlt=pd.DataFrame()
        for name_file in list_file:
            print('Extracting '+ name_file)
            path_file_in=os.path.join(path_in,name_file)
            df_header=ExtractDcmHeader(path_file_in=path_file_in,output_file=False).output
            values=df_header.loc[:,'value']
            #values.reindex([str(num) for num in np.arange(n_keys)])
            values=pd.concat([pd.Series(name_file,index=['file']),values])
            df_header_mlt=df_header_mlt.append(values,ignore_index=True)
        df_header_mlt.to_csv(path_file_out)
        print('Done!')

class ExtractMltDcmHeader():
    def __init__(self,
        path_exp='P:/MRI/pnTTC/Preproc/test_5sub/27_dicom_ses2_t1w/output',
        #path_exp='P:/MRI/pnTTC/Preproc/00_dicom_ses12_exist/pnTTC1_T1',
        #path_exp='P:/MRI/pnTTC/Preproc/00_dicom_ses12_exist/pnTTC2_T1',
        path_file_output='P:/MRI/pnTTC/Preproc/test_5sub/28_header_date_ses2_t1w/output/dcmheader.csv',
        #path_file_output='P:/MRI/pnTTC/Preproc/18_acquisitiondate_t1w/output/dcmheader_ses1.csv',
        #path_file_output='P:/MRI/pnTTC/Preproc/18_acquisitiondate_t1w/output/dcmheader_ses2.csv',
        keys=['Acquisition Date'],
        prefix_dir_sub='CSUB-',
        #suffix_dir_sub='C-01'
        suffix_dir_sub='C-02'
        ):

        list_dir_sub = os.listdir(path_exp)
        list_dir_sub.sort()
        df_header_mlt=pd.DataFrame(columns=['ID_pnTTC','dir_sub','key','value'])
        for dir_sub in list_dir_sub:
            ID_pnTTC=int(dir_sub.replace(prefix_dir_sub,'').replace(suffix_dir_sub,''))
            list_img=os.listdir(os.path.join(path_exp,dir_sub))
            list_img=[img for img in list_img if img.startswith('IM-')]
            list_img.sort()
            path_firstimg=os.path.join(path_exp,dir_sub,list_img[0])
            df_header_firstimg=ExtractDcmHeader(path_firstimg).output
            for key in keys:
                value=df_header_firstimg.loc[df_header_firstimg['name']==key,'value'].values[0]
                df_header_mlt=df_header_mlt.append(pd.Series([ID_pnTTC,dir_sub,key,value],index=df_header_mlt.columns),ignore_index=True)
            print('Checked directory ' + dir_sub + '.')
        
        #df_header_mlt.to_csv(path_file_output,index=False)

        list_id_input=list(df_header_mlt['ID_pnTTC'])
        list_id_output=list(range(1,max(list_id_input)+1,1))
        list_id_diff=list(set(list_id_output)-set(list_id_input))
        list_id_diff.sort()
        df_space=pd.DataFrame(np.nan,columns=df_header_mlt.columns,index=range(len(list_id_diff)))
        df_space['ID_pnTTC']=list_id_diff
        df_header_mlt=pd.concat([df_header_mlt,df_space])
        df_header_mlt=df_header_mlt.sort_values(by='ID_pnTTC')

        df_header_mlt.to_csv(path_file_output,index=False)

        print('All done.')


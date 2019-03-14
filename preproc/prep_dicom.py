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
# make tar.gz file of participant subfolders
##################################################
# used to manipulate pnTTC raw data

class TarGz():
    def __init__(self,
        #path_in='/media/veracrypt2/MRI/pnTTC/Raw/HUMAN-01-ANON_test',
        #path_out='/media/veracrypt1/MRI/pnTTC/Raw/HUMAN-01-ANON_test',
        path_in='/media/veracrypt2/MRI/pnTTC/Raw/HUMAN-01-ANON',
        path_out='/media/veracrypt1/MRI/pnTTC/Raw/HUMAN-01-ANON',
        type_subj=['C-01','C-02','M-01','P-01']
        ):

        print('Starting to pickup subjects and compressing files.')
        list_file = os.listdir(path_in)
        list_file_slctd=[]
        for t_s in type_subj:
            list_file_slctd =list_file_slctd+[f for f in list_file if t_s in f]
        list_file_slctd.sort()
        print('Number of subjects / studies: ' + str(len(list_file_slctd)))

        for f in list_file_slctd:
            path_dir_in=os.path.join(path_in, f)
            path_file_out=os.path.join(path_out,f+'.tar.gz')
            with tarfile.open(path_file_out, "w:gz") as tar:
                tar.add(path_dir_in, arcname=os.path.basename(path_dir_in))
            print('Finished compressiong ' + f)
        
        print('Finished pickup and compressing files.')


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


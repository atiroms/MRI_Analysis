##################################################
# Libraries
##################################################

import os
import shutil
import pandas as pd
import csv
import nilearn.image as nl_image
import json
import numpy as np
import pydicom
import datetime
import gzip


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
# Extract DICOM header metadata
##################################################
# used to extract scan date information

class ExtractDcmHeader():
    def __init__(self,
        path_file_in='P:/MRI/pnTTC/Preproc/00_dicom_ses12_exist/pnTTC2_T1/CSUB-00003C-02/IM-0001-0001-0001.dcm'
        ):

        img_in=pydicom.read_file(path_file_in)
        df_header=pd.DataFrame(columns=['tag','name','value'])
        for element in img_in:
            df_header=df_header.append(pd.Series([str(element.tag),element.name,element.repval],index=df_header.columns),ignore_index=True)
        self.output=df_header
        #print('DICOM element count: '+ str(len(df_header)))
        #print('All done.')

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


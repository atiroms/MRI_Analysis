##################################################
# Libraries
##################################################
import numpy as np
import pandas as pd
import shutil

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
# Extract SPM tissue volume data to combine with CSUB.csv
##################################################
class ExtractTB():
    def __init__(self,
        path_file_src='C:/Users/atiro/Dropbox/MRI_img/pnTTC/puberty/common/tissue_volumes.csv',
        path_file_dst='C:/Users/atiro/Dropbox/MRI_img/pnTTC/puberty/common/tissue_volumes_spaced.csv'
        ):

        # load and calculate global brain calculation data
        df_vol=pd.read_csv(path_file_src,encoding='unicode_escape')
        df_vol['TBV']=df_vol['Volume1']+df_vol['Volume2']
        df_vol['ICV']=df_vol['Volume1']+df_vol['Volume2']+df_vol['Volume3']
        df_vol['ses']=[int(path_file[-16:-14]) for path_file in df_vol.loc[:,'File']]
        df_vol['ID_pnTTC']=[int(path_file[-26:-21]) for path_file in df_vol.loc[:,'File']]

        # split to waves and space
        df_out=pd.DataFrame({'ID_pnTTC':list(range(1,max(df_vol.loc[:,'ID_pnTTC'])+1))})
        df_out['W1_TBV']=pd.Series()
        df_out['W1_ICV']=pd.Series()
        df_out['W2_TBV']=pd.Series()
        df_out['W2_ICV']=pd.Series()
        for idx_row in range(len(df_vol)):
            ses=df_vol.loc[idx_row,'ses']
            id_subj=df_vol.loc[idx_row,'ID_pnTTC']
            if ses==1:
                df_out.loc[df_out['ID_pnTTC']==id_subj,['W1_TBV','W1_ICV']]=df_vol.loc[idx_row,['TBV','ICV']].values.tolist()
            elif ses==2:
                df_out.loc[df_out['ID_pnTTC']==id_subj,['W2_TBV','W2_ICV']]=df_vol.loc[idx_row,['TBV','ICV']].values.tolist()

        # output
        df_out.to_csv(path_file_dst,index=False)

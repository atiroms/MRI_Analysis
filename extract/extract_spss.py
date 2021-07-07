import pandas as pd
import os

path_src='D:/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/common/clinical_data_manipulation'

file_src='181102森田先生.sav'

data_spss=pd.read_spss(os.path.join(path_src,'source',file_src))

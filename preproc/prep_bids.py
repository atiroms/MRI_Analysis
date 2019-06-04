##################################################
# Libraries
##################################################
import os
import shutil
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
# Collect from BIDS folder
##################################################
# collect images from BIDS folder

class Collect():
    def __init__(self,
        #path_src='C:/Users/atiro/Dropbox/temp/BIDS',
        #path_dst='C:/Users/atiro/Dropbox/temp/collect',
        #path_src='C:/Users/NICT_WS/Dropbox/temp/37_c1_bids',
        #path_dst='C:/Users/NICT_WS/Dropbox/temp/54_c1_acpc',
        path_src='C:/Users/NICT_WS/Dropbox/temp/38_c2_bids',
        path_dst='C:/Users/NICT_WS/Dropbox/temp/55_c2_acpc',
        #list_subdir_src=['ses-01/anat']
        list_subdir_src=['ses-02/anat']
        #list_subdir_src=['ses-01/anat','ses-01/fmap']
    ):

        print('Starting Collect()')
        # Create output folder
        print('Starting to create output folder.')
        list_path_mkdir=[]
        list_path_mkdir.append(path_dst)
        list_path_mkdir.append(os.path.join(path_dst,'output'))
        for subdir_src in list_subdir_src:
            list_path_mkdir.append(os.path.join(path_dst,'output',subdir_src))
        for p in list_path_mkdir:
            if not os.path.exists(p):
                os.makedirs(p)
        print('Finished creating output folder.')

        # Copy log folder
        print('Starting to copy log folder.')
        path_log_src=os.path.join(path_src,'log')
        path_log_dst=os.path.join(path_dst,'log')
        shutil.copytree(path_log_src,path_log_dst)
        print('Finished copying log folder.')

        print('Starting to make list of subjects.')
        list_dir_subj = os.listdir(os.path.join(path_src,'output'))
        list_dir_subj = [dir for dir in list_dir_subj if os.path.isdir(os.path.join(path_src,'output',dir))]
        list_dir_subj = [dir for dir in list_dir_subj if dir[:4]=='sub-']        
        list_dir_subj.sort()
        print('Number of subjects: ' + str(len(list_dir_subj)))

        for subdir in list_subdir_src:
            print('Collecting from subdirectory: ' + subdir)
            for dir_subj in list_dir_subj:
                path_subdir=os.path.join(path_src,'output',dir_subj,subdir)
                if os.path.exists(path_subdir):
                    print('Collecting '+ dir_subj)
                    list_file_src=os.listdir(path_subdir)
                    list_file_src=[file_src for file_src in list_file_src if file_src[-7:]=='.nii.gz']
                    list_file_src.sort()
                    for file_src in list_file_src:
                        path_file_src=os.path.join(path_subdir,file_src)
                        file_dst=file_src[:-7]+'.nii'
                        path_file_dst=os.path.join(path_dst,'output',subdir,file_dst)
                        with gzip.open(path_file_src, 'rb') as img_src:
                            with open(path_file_dst, 'wb') as img_dst:
                                shutil.copyfileobj(img_src, img_dst)

        print('Finished Collect()')


##################################################
# Spread back to BIDS folder
##################################################
# spread back images to BIDS folder

class Spread():
    def __init__(self,
        #path_src='C:/Users/atiro/Dropbox/temp/collect',
        #path_dst='C:/Users/atiro/Dropbox/temp/BIDS',
        #path_src='C:/Users/NICT_WS/Dropbox/temp/54_c1_acpc',
        #path_dst='C:/Users/NICT_WS/Dropbox/temp/56_c1_bids',
        path_src='C:/Users/NICT_WS/Dropbox/temp/55_c2_acpc',
        path_dst='C:/Users/NICT_WS/Dropbox/temp/57_c2_bids',
        #list_subdir_src=['ses-01/anat']
        list_subdir_src=['ses-02/anat']
        #list_subdir_src=['ses-01/anat','ses-01/fmap']
    ):

        print('Starting Spread()')
        # Create output folder
        print('Starting to create output folder.')
        list_path_mkdir=[]
        list_path_mkdir.append(path_dst)
        list_path_mkdir.append(os.path.join(path_dst,'output'))
        for p in list_path_mkdir:
            if not os.path.exists(p):
                os.makedirs(p)
        print('Finished creating output folder.')

        # Copy log folder
        #print('Starting to copy log folder.')
        #path_log_src=os.path.join(path_src,'log')
        #path_log_dst=os.path.join(path_dst,'log')
        #shutil.copytree(path_log_src,path_log_dst)
        #print('Finished copying log folder.')

        for subdir in list_subdir_src:
            print('Spreading from subdirectory: ' + subdir)
            list_file_src=os.listdir(os.path.join(path_src,'output',subdir))
            list_file_src.sort()
            for file_src in list_file_src:
                if file_src[-8:]=='_T1w.nii':
                    dir_sub=file_src[:9]
                    dir_ses=file_src[10:16]
                    path_file_src=os.path.join(path_src,'output',subdir,file_src)
                    path_file_dst=os.path.join(path_dst,'output',dir_sub,dir_ses,'anat',file_src+'.gz')
                    os.remove(path_file_dst)
                    with open(path_file_src, 'rb') as img_src:
                        with gzip.open(path_file_dst, 'wb') as img_dst:
                            shutil.copyfileobj(img_src, img_dst)
                    print('Replaced: '+file_src+'.gz')
        print('Finished Spread()')
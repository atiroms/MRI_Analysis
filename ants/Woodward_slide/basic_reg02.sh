#!/bin/bash
source activate neuroimaging
dim=3 # image dimensionality
AP=/home/atiroms/ANTs/bin/ # path to ANTs binaries
ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=4 # controls multi-threading
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
f=/home/atiroms/Dropbox/temp/Shen_to_ICBM/ANTs_practice/data/IXI/T_template2.nii.gz
m=/home/atiroms/Dropbox/temp/Shen_to_ICBM/ANTs_practice/data/IXI594-Guys-1089-T1.nii.gz
mask=/home/atiroms/Dropbox/temp/Shen_to_ICBM/ANTs_practice/data/IXI/T_templateExtractionMask.nii.gz
# fixed and moving image file names and fixed image mask, here the fixed image is the template
if [[ ${#f} -eq 0 ]] ; then #CLI feedback when parameters are not given correctly to the script
  echo usage is
  echo $0 fixed.nii.gz moving.nii.gz fixed_brain_mask.nii.gz
  exit
fi

if [[ ! -s ${nm}0GenericAffine.mat ]] ; then
  $reg -d $dim -r [ $imgs ,1] \
  -m mattes[ $imgs , 1 , 32, regular, 0.05 ] \
  -t translation[ 0.1 ] \
  -c [1000,1.e-8,20] \
  -s 4vox \
  -f 6 -l 1 \
  -m mattes[ $imgs , 1 , 32, regular, 0.1 ] \
  -t rigid[ 0.1 ] \
  -c [1000x1000,1.e-8,20] \
  -s 4x2vox \
  -f 4x2 -l 1 \
  -m mattes[ $imgs , 1 , 32, regular, 0.1 ] \
  -t affine[ 0.1 ] \
  -c [$its,1.e-8,20] \
  -s 4x2x1vox \
  -f 3x2x1 -l 1 \
  -m mattes[ $imgs , 1 , 32 ] \
  -t SyN[ .20, 3, 0 ] \
  -c [ $syn ] \
  -s 1x0.5x0vox \
  -f 4x2x1 -l 1 -u 1 -z 1 -x $mask --float 1 \
  -o [${nm},${nm}_diff.nii.gz,${nm}_inv.nii.gz]
  ${AP}antsApplyTransforms -d $dim -i $m -r $f -n linear -t ${nm}1Warp.nii.gz -t ${nm}0GenericAffine.mat -o ${nm}_warped.nii.gz --float 1
fi

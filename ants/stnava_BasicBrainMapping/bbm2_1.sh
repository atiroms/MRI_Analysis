#!/bin/bash
dim=3 # image dimensionality
AP=/home/atiroms/ANTs/bin/ # /home/yourself/code/ANTS/bin/bin/  # path to ANTs binaries
DATAPATH=/home/atiroms/Dropbox/temp/Shen_to_ICBM/ANTs_practice/data/
OUTPUTPATH=/home/atiroms/Dropbox/temp/Shen_to_ICBM/ANTs_practice/results/
cd ${OUTPUTPATH}
ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=20  # controls multi-threading
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
prefix=synquick_

ImageMath 3 ${DATAPATH}neg_lesion2.nii.gz Neg ${DATAPATH}LesionMap.nii.gz

antsRegistrationSyNQuick.sh \
  -d 3 \
  -m ${DATAPATH}IXI594-Guys-1089-T1.nii.gz \
  -f ${DATAPATH}T1_with_lesion.nii.gz \
  -t s -o ${prefix} \
  -x ${DATAPATH}neg_lesion2.nii.gz

#CreateJacobianDeterminantImage 3 ${nm}_diff1Warp.nii.gz jacobian.nii.gz 1
#CreateWarpedGridImage 3 ${nm}_diff1Warp.nii.gz grid.nii.gz 1x0x1 10x10x10 3x3x3

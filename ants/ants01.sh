#!/bin/bash
source activate neuroimaging
dim=3 # image dimensionality
AP=/home/atiroms/ANTs/bin/ # /home/yourself/code/ANTS/bin/bin/  # path to ANTs binaries
DATAPATH=/home/atiroms/Dropbox/temp/Shen_to_ICBM/image_transform/ants_syn01/source/
OUTPUTPATH=/home/atiroms/Dropbox/temp/Shen_to_ICBM/image_transform/ants_syn01/results2/
cd ${OUTPUTPATH}
ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=20  # controls multi-threading
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
prefix=synquick_

#antsRegistrationSyNQuick.sh \
antsRegistrationSyN.sh \
  -d 3 \
  -m ${DATAPATH}colin27_t1_tal_hires.nii \
  -f ${DATAPATH}mni_icbm152_t1_tal_nlin_asym_09c.nii \
  -t s \
  -o ${prefix} \
  -x ${DATAPATH}mni_icbm152_t1_tal_nlin_asym_09c_mask.nii

antsApplyTransforms \
  -d 3 \
  -i ${DATAPATH}shen_1mm_268_parcellation.nii.gz \
  -r ${DATAPATH}colin27_t1_tal_hires.nii \
  -t ${OUTPUTPATH}synquick_0GenericAffine.mat \
  -t ${OUTPUTPATH}synquick_1Warp.nii.gz \
  -o ${OUTPUTPATH}shen_1mm_268_parcellation_icbm152.nii.gz

sudo docker run -ti --rm -v /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/Nifti3_fMRIPrep:/data:ro -v /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/fMRIPrep_output:/out poldracklab/fmriprep:latest /data /out/out participant --fs-license-file /usr/local/freesurfer/license.txt



##
 fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/Nifti3_fMRIPrep /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/fMRIPrep_output

##
 fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/Nifti3_fMRIPrep /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/fMRIPrep_output
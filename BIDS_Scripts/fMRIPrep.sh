sudo docker run -ti --rm -v /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/Nifti3_fMRIPrep:/data:ro -v /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/fMRIPrep_output:/out poldracklab/fmriprep:latest /data /out/out participant --fs-license-file /usr/local/freesurfer/license.txt



##
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/Nifti3_fMRIPrep /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/fMRIPrep_output


# no tracking and no surface analysis version
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --fs-no-reconall /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/Nifti5_fMRIPrep_nosurf /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/fMRIPrep_nosurf_output



## 
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph /media/veracrypt1/MRI/pnTTC/BIDS/test_2sub/Nifti_1 /media/veracrypt1/MRI/pnTTC/BIDS/test_2sub/fMRIPrep_output_1

fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph /media/veracrypt1/MRI/pnTTC/BIDS/test_2sub/Nifti_2 /media/veracrypt1/MRI/pnTTC/BIDS/test_2sub/fMRIPrep_output_2
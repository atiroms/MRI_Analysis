sudo docker run -ti --rm -v /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/Nifti3_fMRIPrep:/data:ro -v /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/fMRIPrep_output:/out poldracklab/fmriprep:latest /data /out/out participant --fs-license-file /usr/local/freesurfer/license.txt


##
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/Nifti3_fMRIPrep /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/fMRIPrep_output


# no tracking and no surface analysis version
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --fs-no-reconall /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/Nifti5_fMRIPrep_nosurf /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/fMRIPrep_nosurf_output


## 
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph --use-aroma --resource-monitor /media/atiroms/MORITA_HDD4/MRI/pnTTC/BIDS/test_2sub/02_slicetiming /media/atiroms/MORITA_HDD4/MRI/pnTTC/BIDS/test_2sub/03_fmriprep


## 
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph --use-aroma --resource-monitor /media/atiroms/MORITA_HDD4/MRI/pnTTC/BIDS/test_1sub/04_slicetiming_1ses /media/atiroms/MORITA_HDD4/MRI/pnTTC/BIDS/test_1sub/05_fmriprep_1ses

## 
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph --use-aroma --resource-monitor /media/atiroms/MORITA_HDD4/MRI/pnTTC/BIDS/test_1sub/04_slicetiming_1ses /media/atiroms/MORITA_HDD4/MRI/pnTTC/BIDS/test_1sub/05_fmriprep_1ses


## 
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph --use-aroma /media/veracrypt1/MRI/pnTTC/BIDS/03_ses1 /media/veracrypt1/MRI/pnTTC/BIDS/04_fmriprep


## use SyN for susceptibility distortion correction when Fieldmap data is not available
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph --use-aroma --use-syn-sdc /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/04_slicetiming_1ses /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/07_fmriprep_syn_1ses


## no surface analysis, use SyN for susceptibility distortion correction when Fieldmap data is not available
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph --use-aroma --use-syn-sdc --fs-no-reconall /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/04_slicetiming_1ses /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/08_fmriprep_syn_nosurf_1ses


## no surface analysis, for comparison with SyN analysis with no surface analysis
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph --use-aroma --fs-no-reconall /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/04_slicetiming_1ses /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/09_fmriprep_nosurf_1ses


## use SyN for susceptibility distortion correction when Fieldmap data is not available
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph --use-aroma --use-syn-sdc /media/veracrypt1/MRI/pnTTC/BIDS/03_ses1 /media/veracrypt1/MRI/pnTTC/BIDS/05_fmriprep_syn


## use SyN for susceptibility distortion correction when Fieldmap data is not available
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph --use-aroma --use-syn-sdc /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/02_slicetiming_1ses /media/veracrypt1/MRI/pnTTC/BIDS/03_fmriprep_syn


## do not use SyN for susceptibility distortion correction for comparison
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph --use-aroma /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/02_slicetiming_1ses /media/veracrypt1/MRI/pnTTC/BIDS/04_fmriprep_nosyn


## 
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph --use-aroma --bold2t1w-dof=12 --use-syn-sdc /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/04_slicetiming_1ses /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/10_syn_12dof_1ses

## 
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph --use-aroma --use-syn-sdc --force-no-bbr /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/04_slicetiming_1ses /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/11_syn_nobbr_1ses

## 
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph --use-aroma --use-syn-sdc --force-no-bbr --bold2t1w-dof=12 /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/04_slicetiming_1ses /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/12_syn_nobbr_12dof_1ses


## 
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph --use-aroma --bold2t1w-dof=12 --use-syn-sdc /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/02_slicetiming_1ses /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/05_syn_12dof_1ses

## 
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph --use-aroma --use-syn-sdc --force-no-bbr /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/02_slicetiming_1ses /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/06_syn_nobbr_1ses

## 
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph --use-aroma --use-syn-sdc --force-no-bbr --bold2t1w-dof=12 /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/02_slicetiming_1ses /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/07_syn_nobbr_12dof_1ses

## use SyN and 12 dof coregistration
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph --use-aroma --use-syn-sdc --bold2t1w-dof=12 /media/veracrypt1/MRI/pnTTC/BIDS/03_ses1 /media/veracrypt1/MRI/pnTTC/BIDS/07_fmriprep_syn_12dof

## use SyN and 12 dof coregistration, test with 5 sub initial-removed data
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph --use-aroma --use-syn-sdc --bold2t1w-dof=12 /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/09_removeinitial /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/10_remini_syn_12dof

## use SyN and 12 dof coregistration, test with 1 sub initial-removed data, output including T1w native space
fmriprep-docker --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph --use-aroma --use-syn-sdc --output-space T1w template --bold2t1w-dof=12 /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/04_slicetiming_1ses /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/15_syn_12dof_nativeout && echo -e "Subject: Automatic Notification\n\nfMRIPrep on 1 subject done." | sendmail atirom.umusus@gmail.com

## use SyN and 12 dof coregistration, test with 1 sub data, output including T1w native space
fmriprep-docker /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/04_slicetiming_1ses /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/15_syn_12dof_nativeout participant --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid '2mm' --write-graph --use-aroma --use-syn-sdc --output-space T1w template fsaverage5 --bold2t1w-dof=12 && echo -e "Subject: Automatic Notification\n\nfMRIPrep on 1 subject done." | sendmail atirom.umusus@gmail.com

## use SyN and 12 dof coregistration, test with 5 subs data, output including T1w native space
fmriprep-docker /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/02_slicetiming_1ses /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/11_syn_12dof_nativeout participant --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid='2mm' --write-graph --use-aroma --use-syn-sdc --output-space T1w template fsaverage5 --bold2t1w-dof=12 && echo -e "Subject: Automatic Notification\n\nfMRIPrep on 5 subjects done." | sendmail atirom.umusus@gmail.com

## test with 1 sub, use preprocessed 07_recon freesurfer data
fmriprep-docker /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/14_removeinitial /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/21_fmriprep_oldfs participant --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid='1mm' --use-aroma --output-space T1w template fsnative fsaverage --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nfMRIPrep on 1 subject done." | sendmail atirom.umusus@gmail.com

## test with 5 subs, use preprocessed 07_recon freesurfer data
fmriprep-docker /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/09_removeinitial /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/13_fmriprep_oldfs participant --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid='1mm' --use-aroma --output-space T1w template fsnative fsaverage --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nfMRIPrep on 5 subjects done." | sendmail atirom.umusus@gmail.com
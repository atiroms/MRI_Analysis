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
fmriprep-docker /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/09_removeinitial /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/13_fmriprep_oldfs participant --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid='1mm' --use-aroma --output-space T1w template fsnative fsaverage --bold2t1w-dof=6 --nthreads=8 && echo -e "Subject: Automatic Notification\n\nfMRIPrep on 5 subjects done." | sendmail atirom.umusus@gmail.com

## test with 5 subs, use preprocessed 10_recon freesurfer data
fmriprep-docker /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/09_removeinitial /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/16_fmriprep_newfs participant --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid='1mm' --use-aroma --output-space T1w template fsnative fsaverage --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nfMRIPrep on 5 subjects done." | sendmail atirom.umusus@gmail.com

## test with 5 subs, simple procedure
fmriprep-docker /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/09_removeinitial /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/17_fmriprep_nativeout/data participant --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid='1mm' --output-space T1w template fsaverage --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nfMRIPrep on 5 subjects on simple procedure done." | sendmail atirom.umusus@gmail.com

## test with 5 subs, simple procedure re-run
fmriprep-docker /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/09_removeinitial /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/17_fmriprep_simple/output participant --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid='1mm' --output-space T1w template fsaverage --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: test_5sub/17_fmriprep_simple\nStart time: 20190121_2040" | sendmail atirom.umusus@gmail.com

## same as above except ICA-AROMA
fmriprep-docker /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/09_removeinitial /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/20_fmriprep_simple_aroma/output participant --fs-license-file /usr/local/freesurfer/license.txt --notrack --template-resampling-grid='1mm' --use-aroma --output-space T1w template fsaverage --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: test_5sub/20_fmriprep_simple_aroma\nStart time: 20190122_1030" | sendmail atirom.umusus@gmail.com

## re-run from the beginning add --write-graph version 1.2.5
fmriprep-docker /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/24_st_ped/output /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/25_fmriprep/output participant --fs-license-file /usr/local/freesurfer/license.txt --notrack --write-graph --template-resampling-grid='1mm' --use-aroma --output-space T1w template fsaverage --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: test_5sub/25_fmriprep\nStart time: 20190129_1847" | sendmail atirom.umusus@gmail.com

## same as above with singularity image ver 1.2.6-1
# need to copy freesurfer license file to /log
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/26_fmriprep_latest:${HOME}/data /data/applications/fmriprep-1261.simg ${HOME}/data/input ${HOME}/data/output participant --fs-license-file ${HOME}/data/log/license.txt --notrack --write-graph --template-resampling-grid='1mm' --use-aroma --output-space T1w template fsaverage --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: test_5sub/26_fmriprep_latest\nStart time: 20190129_2108" | sendmail atirom.umusus@gmail.com

## same as above except using SyN
# need to copy freesurfer license file to /log
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/29_fmriprep_latest_syn:${HOME}/data /data/applications/fmriprep-1261.simg ${HOME}/data/input ${HOME}/data/output participant --fs-license-file ${HOME}/data/log/license.txt --notrack --write-graph --template-resampling-grid='1mm' --use-syn-sdc --use-aroma --output-space T1w template fsaverage --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: test_5sub/29_fmriprep_latest_syn\nStart time: 20190131_0841" | sendmail atirom.umusus@gmail.com

## same as above except output without T1w space
# need to copy freesurfer license file to /log
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/31_fmriprep_latest_syn_templateout:${HOME}/data /data/applications/fmriprep-1261.simg ${HOME}/data/input ${HOME}/data/output participant --fs-license-file ${HOME}/data/log/license.txt --notrack --write-graph --template-resampling-grid='1mm' --use-syn-sdc --use-aroma --output-space template fsaverage --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: test_5sub/31_fmriprep_latest_syn_templateout\nStart time: 20190131_1806" | sendmail atirom.umusus@gmail.com

## same as above except specifying working directory and output only in native space
# need to copy freesurfer license file to /log
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/34_fmriprep_latest_syn_nativeout:${HOME}/data /data/applications/fmriprep-1261.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --write-graph --template-resampling-grid='1mm' --use-syn-sdc --use-aroma --output-space T1w --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: test_5sub/34_fmriprep_latest_syn_nativeout\nStart time: 20190204_1630" | sendmail atirom.umusus@gmail.com

## same as above except output resampling grid in 2mm and output only in template space
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/35_fmriprep_latest_syn_templateout_2mm:${HOME}/data /data/applications/fmriprep-1261.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --write-graph --template-resampling-grid='2mm' --use-syn-sdc --use-aroma --output-space template --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: test_5sub/35_fmriprep_latest_syn_templateout_2mm\nStart time: 20190204_2030" | sendmail atirom.umusus@gmail.com

## same as above except output resampling grid in 2mm and output only in template space
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/36_fmriprep_latest_syn_templateout_native:${HOME}/data /data/applications/fmriprep-1261.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --write-graph --template-resampling-grid='native' --use-syn-sdc --use-aroma --output-space template --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: test_5sub/36_fmriprep_latest_syn_templateout_native\nStart time: 20190204_2030" | sendmail atirom.umusus@gmail.com

## run on real data (at last!)
# first half
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/19_1_fmriprep:${HOME}/data /data/applications/fmriprep-1261.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --use-syn-sdc --use-aroma --output-space T1w template --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 19_1_fmriprep\nStart time: 20190212_2245" | sendmail atirom.umusus@gmail.com

# second half
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/19_2_fmriprep:${HOME}/data /data/applications/fmriprep-1261.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --use-syn-sdc --use-aroma --output-space T1w template --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 19_2_fmriprep\nStart time: 20190212_2245" | sendmail atirom.umusus@gmail.com


## test with fieldmap for 1 sub
# 42 use fieldmap, 6 dof, without syn
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_1sub/42_fmriprep_fmap:${HOME}/data /data/applications/fmriprep-132.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --use-aroma --output-space T1w template --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 42_fmriprep_fmap\nStart time: 20190403_1730" | sendmail atirom.umusus@gmail.com

# 43 ignore fieldmap, 6 dof, without syn
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_1sub/43_fmriprep:${HOME}/data /data/applications/fmriprep-132.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --ignore fieldmaps --use-aroma --output-space T1w template --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 43_fmriprep\nStart time: 20190403_1730" | sendmail atirom.umusus@gmail.com

# 44 use fieldmap, 6 dof, force syn
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_1sub/44_fmriprep_fmap_syn:${HOME}/data /data/applications/fmriprep-132.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --force-syn --use-aroma --output-space T1w template --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 44_fmriprep_fmap_syn\nStart time: 20190403_1730" | sendmail atirom.umusus@gmail.com

# 45 ignore fieldmap, 6 dof, use syn
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_1sub/45_fmriprep_syn:${HOME}/data /data/applications/fmriprep-132.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --ignore fieldmaps --use-syn-sdc --use-aroma --output-space T1w template --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 45_fmriprep_syn\nStart time: 20190403_1730" | sendmail atirom.umusus@gmail.com

# 46 use fieldmap, 9 dof, force syn
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_1sub/46_fmriprep_fmap_syn_9dof:${HOME}/data /data/applications/fmriprep-132.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --force-syn --use-aroma --output-space T1w template --bold2t1w-dof=9 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 46_fmriprep_fmap_syn_9dof\nStart time: 20190403_2000" | sendmail atirom.umusus@gmail.com

# 47 ignore fieldmap, 9 dof, use syn
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_1sub/47_fmriprep_syn_9dof:${HOME}/data /data/applications/fmriprep-132.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --ignore fieldmaps --use-syn-sdc --use-aroma --output-space T1w template --bold2t1w-dof=9 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 47_fmriprep_syn_9dof\nStart time: 20190403_2000" | sendmail atirom.umusus@gmail.com

# 48 use fieldmap, 12 dof, force syn
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_1sub/48_fmriprep_fmap_syn_12dof:${HOME}/data /data/applications/fmriprep-132.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --force-syn --use-aroma --output-space T1w template --bold2t1w-dof=12 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 48_fmriprep_fmap_syn_12dof\nStart time: 20190403_2000" | sendmail atirom.umusus@gmail.com

# 49 ignore fieldmap, 12 dof, use syn
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_1sub/49_fmriprep_syn_12dof:${HOME}/data /data/applications/fmriprep-132.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --ignore fieldmaps --use-syn-sdc --use-aroma --output-space T1w template --bold2t1w-dof=12 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 49_fmriprep_syn_12dof\nStart time: 20190403_2000" | sendmail atirom.umusus@gmail.com


## test with 5 subjects, revert to fmriprep 1.2.6-1
# 55_01 ignore fieldmap, 6 dof, no SyN
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/55_01_fmriprep:${HOME}/data /data/applications/fmriprep-1261.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --ignore fieldmaps --use-aroma --output-space T1w template --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 55_01_fmriprep\nStart time: 20190404_1330" | sendmail atirom.umusus@gmail.com

# 55_02 use fieldmap, 6 dof, no SyN
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/55_02_fmriprep:${HOME}/data /data/applications/fmriprep-1261.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --use-aroma --output-space T1w template --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 55_02_fmriprep\nStart time: 20190404_1330" | sendmail atirom.umusus@gmail.com

# 55_03 ignore fieldmap, 6 dof, use SyN
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/55_03_fmriprep:${HOME}/data /data/applications/fmriprep-1261.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --ignore fieldmaps --use-syn-sdc --use-aroma --output-space T1w template --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 55_03_fmriprep\nStart time: 20190404_1330" | sendmail atirom.umusus@gmail.com

# 55_04 use fieldmap, 6 dof, force SyN
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/55_04_fmriprep:${HOME}/data /data/applications/fmriprep-1261.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --force-syn --use-aroma --output-space T1w template --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 55_04_fmriprep\nStart time: 20190404_1330" | sendmail atirom.umusus@gmail.com

# 55_05 ignore fieldmap, 9 dof, no SyN
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/55_05_fmriprep:${HOME}/data /data/applications/fmriprep-1261.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --ignore fieldmaps --use-aroma --output-space T1w template --bold2t1w-dof=9 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 55_05_fmriprep\nStart time: 20190404_1330" | sendmail atirom.umusus@gmail.com

# 55_06 use fieldmap, 9 dof, no SyN
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/55_06_fmriprep:${HOME}/data /data/applications/fmriprep-1261.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --use-aroma --output-space T1w template --bold2t1w-dof=9 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 55_06_fmriprep\nStart time: 20190404_1330" | sendmail atirom.umusus@gmail.com

# 55_07 ignore fieldmap, 9 dof, use SyN
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/55_07_fmriprep:${HOME}/data /data/applications/fmriprep-1261.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --ignore fieldmaps --use-syn-sdc --use-aroma --output-space T1w template --bold2t1w-dof=9 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 55_07_fmriprep\nStart time: 20190404_1330" | sendmail atirom.umusus@gmail.com

# 55_08 use fieldmap, 9 dof, force SyN
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/55_08_fmriprep:${HOME}/data /data/applications/fmriprep-1261.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --force-syn --use-aroma --output-space T1w template --bold2t1w-dof=9 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 55_08_fmriprep\nStart time: 20190404_1330" | sendmail atirom.umusus@gmail.com

# 55_09 ignore fieldmap, 12 dof, no SyN
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/55_09_fmriprep:${HOME}/data /data/applications/fmriprep-1261.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --ignore fieldmaps --use-aroma --output-space T1w template --bold2t1w-dof=12 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 55_09_fmriprep\nStart time: 20190404_1330" | sendmail atirom.umusus@gmail.com

# 55_10 use fieldmap, 12 dof, no SyN
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/55_10_fmriprep:${HOME}/data /data/applications/fmriprep-1261.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --use-aroma --output-space T1w template --bold2t1w-dof=12 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 55_10_fmriprep\nStart time: 20190404_1330" | sendmail atirom.umusus@gmail.com

# 55_11 ignore fieldmap, 12 dof, use SyN
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/55_11_fmriprep:${HOME}/data /data/applications/fmriprep-1261.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --ignore fieldmaps --use-syn-sdc --use-aroma --output-space T1w template --bold2t1w-dof=12 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 55_11_fmriprep\nStart time: 20190404_1330" | sendmail atirom.umusus@gmail.com

# 55_12 use fieldmap, 12 dof, force SyN
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/55_12_fmriprep:${HOME}/data /data/applications/fmriprep-1261.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --force-syn --use-aroma --output-space T1w template --bold2t1w-dof=12 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 55_12_fmriprep\nStart time: 20190404_1330" | sendmail atirom.umusus@gmail.com


## analysis on real data
## first wave
# 41_c1_1
singularity run --cleanenv -B /media/veracrypt2/MRI/pnTTC/Preproc/41_c1_1_fmriprep:${HOME}/data /data/applications/fmriprep-132.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --use-aroma --output-space T1w template --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 41_c1_1_fmriprep\nStart time: 20190418_1700" | sendmail atirom.umusus@gmail.com

# 41_c1_2
singularity run --cleanenv -B /media/veracrypt2/MRI/pnTTC/Preproc/41_c1_2_fmriprep:${HOME}/data /data/applications/fmriprep-132.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --use-aroma --output-space T1w template --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 41_c1_2_fmriprep\nStart time: 20190424_0900" | sendmail atirom.umusus@gmail.com

## second wave
# 42_c2_1
singularity run --cleanenv -B /media/veracrypt2/MRI/pnTTC/Preproc/42_c2_1_fmriprep:${HOME}/data /data/applications/fmriprep-132.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --use-aroma --output-space T1w template --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 42_c2_1_fmriprep\nStart time: 201904**_****" | sendmail atirom.umusus@gmail.com

# 42_c2_2
singularity run --cleanenv -B /media/veracrypt2/MRI/pnTTC/Preproc/42_c2_2_fmriprep:${HOME}/data /data/applications/fmriprep-132.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --use-aroma --output-space T1w template --bold2t1w-dof=6 && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: 42_c2_2_fmriprep\nStart time: 20190424_0900" | sendmail atirom.umusus@gmail.com

## re-run on ac-pc coregistered data
# c1
singularity run --cleanenv -B /media/veracrypt2/MRI_img/pnTTC/preproc/60_c1_fmriprep:${HOME}/data /data/applications/fmriprep-140.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --use-aroma --output-space T1w template --bold2t1w-dof=6

# c2
singularity run --cleanenv -B /media/veracrypt2/MRI_img/pnTTC/preproc/61_c2_fmriprep:${HOME}/data /data/applications/fmriprep-140.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --use-aroma --output-space T1w template --bold2t1w-dof=6

## re-run on ac-pc coregistered data, ignore fieldmaps
# fmriprep-132  on Ubuntu_2
singularity run --cleanenv -B /media/veracrypt2/MRI_img/pnTTC/preproc/test_fmriprep_1:${HOME}/data /data/applications/fmriprep-132.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --ignore fieldmaps --use-syn-sdc --use-aroma --output-space T1w template --bold2t1w-dof=6

# fmriprep-140 on Ubuntu_1
singularity run --cleanenv -B /media/veracrypt1/MRI_img/pnTTC/preproc/test_fmriprep_2:${HOME}/data /data/applications/fmriprep-140.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --ignore fieldmaps --use-syn-sdc --use-aroma --output-space T1w template --bold2t1w-dof=6

## test on 3 subjects
# 01
singularity run --cleanenv -B /media/veracrypt1/MRI_img/pnTTC/test3/01_fmriprep:${HOME}/data /data/applications/fmriprep-132.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --ignore fieldmaps --use-syn-sdc --use-aroma --output-space T1w template --bold2t1w-dof=6

# 02
singularity run --cleanenv -B /media/veracrypt1/MRI_img/pnTTC/test3/02_fmriprep:${HOME}/data /data/applications/fmriprep-132.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --ignore fieldmaps --use-syn-sdc --use-aroma --output-space T1w template --bold2t1w-dof=6

# 03
singularity run --cleanenv -B /media/veracrypt1/MRI_img/pnTTC/test3/03_fmriprep:${HOME}/data /data/applications/fmriprep-140.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --ignore fieldmaps --use-syn-sdc --use-aroma --output-space T1w template --bold2t1w-dof=6

# 04
singularity run --cleanenv -B /media/veracrypt1/MRI_img/pnTTC/test3/04_fmriprep:${HOME}/data /data/applications/fmriprep-140.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --ignore fieldmaps --use-syn-sdc --use-aroma --output-space T1w template --bold2t1w-dof=6

## 67_c1_fmriprep
singularity run --cleanenv -B /media/veracrypt2/MRI_img/pnTTC/preproc/67_c1_fmriprep:${HOME}/data /data/applications/fmriprep-132.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --ignore fieldmaps --use-syn-sdc --use-aroma --output-space T1w template --bold2t1w-dof=6

## 68_c2_fmriprep
singularity run --cleanenv -B /media/veracrypt2/MRI_img/pnTTC/preproc/68_c2_fmriprep:${HOME}/data /data/applications/fmriprep-132.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --template-resampling-grid='2mm' --ignore fieldmaps --use-syn-sdc --use-aroma --output-space T1w template --bold2t1w-dof=6

# 03_5sub_c1_fmriprep
singularity run --cleanenv -B /media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/5sub/03_5sub_c1_fmriprep:${HOME}/data /data/applications/fmriprep-154.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --ignore fieldmaps --use-syn-sdc --use-aroma --output-spaces MNI152NLin2009cAsym:res-2 MNI152NLin6Asym:res-2 T1w --bold2t1w-dof=6

# 04_5sub_c1_fmriprep
singularity run --cleanenv -B /media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/5sub/04_5sub_c1_fmriprep:${HOME}/data /data/applications/fmriprep-154.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --ignore fieldmaps --use-syn-sdc --use-aroma --output-spaces MNI152NLin2009cAsym:res-2 --bold2t1w-dof=6

## 401_c1_fmriprep
singularity run --cleanenv -B /media/veracrypt2/MRI_img/pnTTC/preproc/401_c1_fmriprep:${HOME}/data /data/applications/fmriprep-154.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --ignore fieldmaps --use-syn-sdc --use-aroma --output-spaces MNI152NLin2009cAsym:res-2 MNI152NLin6Asym:res-2 --bold2t1w-dof=6

## 402_c2_fmriprep
singularity run --cleanenv -B /media/veracrypt2/MRI_img/pnTTC/preproc/402_c2_fmriprep:${HOME}/data /data/applications/fmriprep-154.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --ignore fieldmaps --use-syn-sdc --use-aroma --output-spaces MNI152NLin2009cAsym:res-2 MNI152NLin6Asym:res-2 --bold2t1w-dof=6

## 05_5sub_c1_fmriprep
singularity run --cleanenv -B /media/atiroms/SSD_2TB/MRI_img/pnTTC/preproc/5sub/05_5sub_c1_fmriprep:${HOME}/data /data/applications/fmriprep-154.simg ${HOME}/data/input ${HOME}/data/output participant --work-dir ${HOME}/data/output/work --fs-license-file ${HOME}/data/log/license.txt --notrack --ignore fieldmaps --use-syn-sdc --use-aroma --output-spaces MNI152NLin2009cAsym:res-2 MNI152NLin6Asym:res-2 --bold2t1w-dof=6
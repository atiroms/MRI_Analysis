sudo docker run --rm -it -v /media/atiroms/MORITA_HDD4/MRI/pnTTC/BIDS/test:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC1_T1/CSUB-{subject}C-01/*.dcm -o /base/Nifti/ -f convertall -s 00003 00004 -c none --overwrite


sudo docker run --rm -it -v /media/atiroms/MORITA_HDD4/MRI/pnTTC/BIDS/test:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC1_*/CSUB-{subject}C-01/*.dcm -o /base/Nifti/ -f convertall -s 00003 00004 -c none --overwrite


sudo docker run --rm -it -v /media/veracrypt1/MRI/pnTTC/BIDS/test2:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC1_T1/CSUB-{subject}C-01/*.dcm -o /base/Nifti/ -f convertall -s 00014 00019 -ss 01 -c none --overwrite


sudo docker run --rm -it -v /media/veracrypt1/MRI/pnTTC/BIDS/test2:/base nipy/heudiconv:latest -d /base/Dicom/pnTTC1_rsfMRI/CSUB-{subject}C-01/*.dcm -o /base/Nifti/ -f convertall -s 00014 00019 -ss 01 -c none --overwrite


# test on raw data
sudo docker run --rm -it -v /media/veracrypt1/MRI/pnTTC/Preproc/test_1sub/30_heudiconv:/base nipy/heudiconv:latest -d /base/input/CSUB-{subject}C-01/*/*/*/*.dcm -o /base/output/ -f convertall -s 00014 -ss 01 -c none --overwrite
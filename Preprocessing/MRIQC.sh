##
sudo docker run -it --rm -v /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/Nifti4_MRIQC:/data:ro -v /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/MRIQC_output:/out poldracklab/mriqc:latest /data /out participant --no-sub


##
sudo docker run -it --rm -v /media/atiroms/MORITA_HDD4/MRI/pnTTC/BIDS/test_2sub/02_slicetiming:/data:ro -v /media/atiroms/MORITA_HDD4/MRI/pnTTC/BIDS/test_2sub/04_mriqc:/out poldracklab/mriqc:latest /data /out participant --no-sub


##
sudo docker run -it --rm -v /media/atiroms/MORITA_HDD4/MRI/pnTTC/BIDS/test_1sub/04_slicetiming_1ses:/data:ro -v /media/atiroms/MORITA_HDD4/MRI/pnTTC/BIDS/test_1sub/06_mriqc_1ses:/out poldracklab/mriqc:latest /data /out participant --no-sub


##
source activate neuroimaging
sudo docker run -it --rm -v /media/veracrypt1/MRI/pnTTC/BIDS/03_ses1:/data:ro -v /media/veracrypt1/MRI/pnTTC/BIDS/05_mriqc:/out poldracklab/mriqc:latest /data /out participant --no-sub

# MRIQC on all subjects with T1 existance, with verbose results
# ses1
source activate neuroimaging
sudo docker run -it --rm -v /media/veracrypt1/MRI/pnTTC/BIDS/10_bids_ses1_t1exist:/data:ro -v /media/veracrypt1/MRI/pnTTC/BIDS/12_mriqc_ses1_t1exist:/out poldracklab/mriqc:latest /data /out participant --no-sub --verbose-reports && echo -e "Subject: Automatic Notification\n\nMRIQC on all ses1 subjects with T1 images done." | sendmail atirom.umusus@gmail.com

# ses2
source activate neuroimaging
sudo docker run -it --rm -v /media/veracrypt1/MRI/pnTTC/BIDS/11_bids_ses2_t1exist:/data:ro -v /media/veracrypt1/MRI/pnTTC/BIDS/13_mriqc_ses2_t1exist:/out poldracklab/mriqc:latest /data /out participant --no-sub --verbose-reports && echo -e "Subject: Automatic Notification\n\nMRIQC on all ses2 subjects with T1 images done." | sendmail atirom.umusus@gmail.com
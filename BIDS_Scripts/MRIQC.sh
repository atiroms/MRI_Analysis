##
sudo docker run -it --rm -v /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/Nifti4_MRIQC:/data:ro -v /media/veracrypt1/MRI/pnTTC/BIDS/test_3sub/MRIQC_output:/out poldracklab/mriqc:latest /data /out participant --no-sub


##
sudo docker run -it --rm -v /media/atiroms/MORITA_HDD4/MRI/pnTTC/BIDS/test_2sub/02_slicetiming:/data:ro -v /media/atiroms/MORITA_HDD4/MRI/pnTTC/BIDS/test_2sub/04_mriqc:/out poldracklab/mriqc:latest /data /out participant --no-sub


##
sudo docker run -it --rm -v /media/atiroms/MORITA_HDD4/MRI/pnTTC/BIDS/test_1sub/04_slicetiming_1ses:/data:ro -v /media/atiroms/MORITA_HDD4/MRI/pnTTC/BIDS/test_1sub/06_mriqc_1ses:/out poldracklab/mriqc:latest /data /out participant --no-sub
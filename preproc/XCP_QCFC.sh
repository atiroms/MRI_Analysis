## test running QCFC of XCP after manually integrating XCP parallel-processed results
singularity run --cleanenv -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/45_xcp_qcfc:${HOME}/data /data/applications/xcpEngine-070-20190130.simg -d ${HOME}/data/input/fc-acompcor_fconly_noqcfc_power.dsn -c ${HOME}/data/input/func_cohort.csv -o ${HOME}/data/output -t 1 -r ${HOME}/data && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: test_5sub/45_xcp_qcfc\nStart time: 20190218_1120" | sendmail atirom.umusus@gmail.com
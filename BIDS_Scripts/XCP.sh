## test anatomical analysis on 1 subject
singularity run -B /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/15_xcp_anat:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/anat-antsct.dsn -c ${HOME}/data/anat_cohort.csv -o ${HOME}/data/xcp_output -t 1 -r ${HOME}/data

## test functional analysis on 1 subject
singularity run -B /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/16_xcp_func:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/fc-36p.dsn -c ${HOME}/data/func_cohort.csv -o ${HOME}/data/xcp_output -t 1 -r ${HOME}/data

## test anatomical and functional analysis on 1 subject and send email
singularity run -B /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/15_xcp:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/anat-antsct.dsn -c ${HOME}/data/anat_cohort.csv -o ${HOME}/data/xcp_output -t 1 -r ${HOME}/data && singularity run -B /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/15_xcp:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/fc-36p.dsn -c ${HOME}/data/func_cohort.csv -o ${HOME}/data/xcp_output -t 1 -r ${HOME}/data && echo -e "Subject: Automatic Notification\n\nAnalysis done" | sendmail atirom.umusus@gmail.com
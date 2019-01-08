## test anatomical analysis on 1 subject
singularity run -B /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/15_xcp_anat:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/anat-antsct.dsn -c ${HOME}/data/anat_cohort.csv -o ${HOME}/data/xcp_output -t 1 -r ${HOME}/data

## test functional analysis on 1 subject
singularity run -B /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/16_xcp_func:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/fc-36p.dsn -c ${HOME}/data/func_cohort.csv -o ${HOME}/data/xcp_output -t 1 -r ${HOME}/data

## test anatomical and functional analysis on 1 subject and send email
singularity run -B /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/15_xcp:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/anat-antsct.dsn -c ${HOME}/data/anat_cohort.csv -o ${HOME}/data/xcp_output -t 1 -r ${HOME}/data && singularity run -B /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/15_xcp:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/fc-36p.dsn -c ${HOME}/data/func_cohort.csv -o ${HOME}/data/xcp_output -t 1 -r ${HOME}/data && echo -e "Subject: Automatic Notification\n\nAnat and Func analysis on 1 subject done" | sendmail atirom.umusus@gmail.com

## test anatomical and functional analysis on 1 subject and send email (tutorial dataset)
singularity run -B /media/veracrypt1/MRI/XCP_tutorial:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/anat-antsct.dsn -c ${HOME}/data/anat_cohort.csv -o ${HOME}/data/xcp_output -t 1 -r ${HOME}/data && singularity run -B /media/veracrypt1/MRI/XCP_tutorial:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/fc-36p.dsn -c ${HOME}/data/func_cohort.csv -o ${HOME}/data/xcp_output -t 1 -r ${HOME}/data && echo -e "Subject: Automatic Notification\n\nAnat and Func analysis on tutorial dataset done" | sendmail atirom.umusus@gmail.com

## functional-only analysis on 1 subject and send email (original dataset)
singularity run -B /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/16_xcp:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/fc-36p_dvo.dsn -c ${HOME}/data/func_cohort.csv -o ${HOME}/data/xcp_output -t 1 -r ${HOME}/data && echo -e "Subject: Automatic Notification\n\nFunc-only analysis on original dataset of 1 sub done" | sendmail atirom.umusus@gmail.com

## functional-only analysis on 1 subject and send email (original dataset) (re-run with /anat contents copied to session-specific directory) * error with unmatched volume number, probably caused by volume deletion
singularity run -B /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/17_xcp:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/fc-36p_dvo.dsn -c ${HOME}/data/func_cohort.csv -o ${HOME}/data/xcp_output -t 1 -r ${HOME}/data && echo -e "Subject: Automatic Notification\n\nFunc-only analysis on original dataset of 1 sub done" | sendmail atirom.umusus@gmail.com

## functional-only analysis on 1 subject and send email (original dataset) (re-run with /anat contents copied to session-specific directory)
singularity run -B /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/18_xcp:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/fc-36p.dsn -c ${HOME}/data/func_cohort.csv -o ${HOME}/data/xcp_output -t 1 -r ${HOME}/data && echo -e "Subject: Automatic Notification\n\nFunc-only analysis on original dataset of 1 sub done" | sendmail atirom.umusus@gmail.com

## functional-only analysis on 1 subject and send email (original dataset) (re-run with /anat contents copied to session-specific directory) * only coregistration metrics copied *error
singularity run -B /media/veracrypt1/MRI/pnTTC/BIDS/test_1sub/19_xcp:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/fc-36p.dsn -c ${HOME}/data/func_cohort.csv -o ${HOME}/data/xcp_output -t 1 -r ${HOME}/data && echo -e "Subject: Automatic Notification\n\nFunc-only analysis on original dataset of 1 sub done" | sendmail atirom.umusus@gmail.com
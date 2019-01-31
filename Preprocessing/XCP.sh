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

## functional-only analysis on 5 subjects and send email (original dataset) (/anat contents moved to session-specific directory using MoveAnat class)
singularity run -B /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/12_xcp:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/fc-36p.dsn -c ${HOME}/data/func_cohort.csv -o ${HOME}/data/xcp_output -t 1 -r ${HOME}/data && echo -e "Subject: Automatic Notification\n\nFunc-only analysis on original dataset of 5 subs done." | sendmail atirom.umusus@gmail.com

## func-only analysis, FC only output
singularity run -B /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/15_xcp_fconly:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/fc-36p_fconly.dsn -c ${HOME}/data/func_cohort.csv -o ${HOME}/data/xcp_output -t 1 -r ${HOME}/data && echo -e "Subject: Automatic Notification\n\nXCP func-only analysis on original dataset of 5 subs (functional connectivity only) done." | sendmail atirom.umusus@gmail.com

## comparison of template input and native input
# template input
singularity run -B /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/18_xcp_templatein/analysis:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/fc-36p_fconly.dsn -c ${HOME}/data/func_cohort.csv -o ${HOME}/data/xcp_output -t 1 -r ${HOME}/data && echo -e "Subject: Automatic Notification\n\nXCP analysis using template input done (for comparison with native input)." | sendmail atirom.umusus@gmail.com

# native input
singularity run -B /media/veracrypt1/MRI/pnTTC/BIDS/test_5sub/19_xcp_nativein/analysis:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/fc-36p_fconly.dsn -c ${HOME}/data/func_cohort.csv -o ${HOME}/data/xcp_output -t 1 -r ${HOME}/data && echo -e "Subject: Automatic Notification\n\nXCP analysis using native input done (for comparison with template input)." | sendmail atirom.umusus@gmail.com

# template input re-run
singularity run -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/18_xcp_templatein:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/input/fc-36p_fconly.dsn -c ${HOME}/data/input/func_cohort.csv -o ${HOME}/data/output -t 1 -r ${HOME}/data && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: test_5sub/18_xcp_templatein\nStart time: 20190121_2224" | sendmail atirom.umusus@gmail.com

# native input re-run
singularity run -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/19_xcp_nativein:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/input/fc-36p_fconly.dsn -c ${HOME}/data/input/func_cohort.csv -o ${HOME}/data/output -t 1 -r ${HOME}/data && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: test_5sub/19_xcp_nativein\nStart time: 20190121_2224" | sendmail atirom.umusus@gmail.com

## comparison of aroma output with and without fMRIPrep --use-aroma 
# without fMRIPrep ICA-AROMA
singularity run -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/21_xcp_aroma:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/input/fc-aroma_fconly.dsn -c ${HOME}/data/input/func_cohort.csv -o ${HOME}/data/output -t 1 -r ${HOME}/data && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: test_5sub/21_xcp_aroma\nStart time: 20190122_1850" | sendmail atirom.umusus@gmail.com

# with fMRIPrep ICA-AROMA
singularity run -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/22_xcp_aroma_aromain:${HOME}/data /data/applications/xcpEngine.simg -d ${HOME}/data/input/fc-aroma_fconly.dsn -c ${HOME}/data/input/func_cohort.csv -o ${HOME}/data/output -t 1 -r ${HOME}/data && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: test_5sub/22_xcp_aroma_aromain\nStart time: 20190122_1850" | sendmail atirom.umusus@gmail.com

## re-run with new fmriprep data, latest singularity image
singularity run -B /media/veracrypt1/MRI/pnTTC/Preproc/test_5sub/30_xcp_36p:${HOME}/data /data/applications/xcpEngine-070-20190130.simg -d ${HOME}/data/input/fc-36p_fconly.dsn -c ${HOME}/data/input/func_cohort.csv -o ${HOME}/data/output -t 1 -r ${HOME}/data && echo -e "Subject: Automatic Notification\n\nAutomatic notification of analysis completion.\n\nAnalysis: test_5sub/30_xcp_36p\nStart time: 20190131_1714" | sendmail atirom.umusus@gmail.com
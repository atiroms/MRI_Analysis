% List of open inputs
% Factorial design specification: Directory - cfg_files
% Factorial design specification: Scans - cfg_files
% Factorial design specification: Scans - cfg_files
% Factorial design specification: Scans - cfg_files
% Factorial design specification: Scans - cfg_files
%nrun = 8; % enter the number of runs here


%path_home='P:\MRI\pnTTC\Prosociality_DC_Dr_Okada\SPM\SPM2';
%path_clinicaldata='P:\MRI\pnTTC\Prosociality_DC_Dr_Okada\Info\ClinicalData.csv';
%jobfile = {'D:\atiroms\GitHub\MRI_Analysis\SPM\batch2_anova\batch01_job.m'};

%path_home='/media/veracrypt2/MRI/pnTTC/Prosociality_DC_Dr_Okada/SPM/SPM2';
%path_home='/media/veracrypt2/MRI/pnTTC/Prosociality_DC_Dr_Okada/SPM/SPM3';
path_home='/media/veracrypt2/MRI/pnTTC/Prosociality_DC_Dr_Okada/SPM/SPM5';
%path_clinicaldata='/media/veracrypt2/MRI/pnTTC/Prosociality_DC_Dr_Okada/Info/ClinicalData.csv';
path_clinicaldata='/media/veracrypt2/MRI/pnTTC/Prosociality_DC_Dr_Okada/Info/ClinicalData.csv';
clinical_data=readtable(path_clinicaldata);

rois=[1,2,3,4,5,6,7,8];
n_rois=length(rois);

models=[1,2,3];
%model=2;

sibling_groups={'1.Only', '2.First-born', '3.Middle-born', '4.Last-born'};
n_sibling_groups=length(sibling_groups);


for model=models
    %jobfile = {'/home/atiroms/Documents/GitHub/MRI_Analysis/SPM/batch2_anova/batch02_model2_job.m'};
    jobfile = {join(['/home/atiroms/Documents/GitHub/MRI_Analysis/SPM/batch2_anova/batch02_model',num2str(model),'_job.m'])};
    jobs = repmat(jobfile, 1, n_rois);
    inputs = cell(5, n_rois);
    
    for cnt_roi = 1:n_rois
        roi=rois(cnt_roi);
    
        path_out=fullfile(path_home,join(['Model',num2str(model),'_ROI',num2str(roi)]));
        if exist(path_out)~=7
            mkdir(path_out);
        end
        inputs{1, cnt_roi} = {path_out}; % Factorial design specification: Directory - cfg_files
        
        for cnt_sibling_group=1:n_sibling_groups
            scans={};
            cnt_data=0;
            for i=1:height(clinical_data)
                if clinical_data{i,join(['Model',num2str(model)])}==1
                    if clinical_data{i,'Sibling_status'}{1}==string(sibling_groups(cnt_sibling_group))
                        cnt_data=cnt_data+1;
                        id_pnttc=num2str(clinical_data{i,'ID'},'%5.5u');
                        path_zroifc=fullfile(path_home,'FC_FunImgARglobalCWSF',join(['zROI',num2str(roi),'FCMap_CSUB-',id_pnttc,'C-01.nii,1']));
                        scans(cnt_data,1)={path_zroifc};
                    end
                end
            end
            inputs{cnt_sibling_group+1,cnt_roi}=scans;
        end
    end
    
    spm('defaults', 'PET');
    spm_jobman('run', jobs, inputs{:});
    
end
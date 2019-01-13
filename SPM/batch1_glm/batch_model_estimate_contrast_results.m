% List of open inputs
% Factorial design specification: Directory - cfg_files
% Factorial design specification: Scans - cfg_files

%path_home='H:\MRI\pnTTC\Prosociality_DC_Dr_Okada\SPM';
path_home='/media/veracrypt2/MRI/pnTTC/Prosociality_DC_Dr_Okada/SPM/SPM1';

rois=[1,2,3,4,5,6,7,8];
models=[1,2,3];
n_rois=length(rois);

for model=models
    clinical_data=csvread(fullfile(path_home,'Script',join(['model',num2str(model),'.csv'])),1,0);

    jobfile = {fullfile(path_home,'Script',join(['batch_model_estimate_contrast_results_model',num2str(model),'_job.m']))};

    jobs = repmat(jobfile, 1, n_rois);
    
    inputs = cell(2, n_rois);

    for cnt_roi = 1:n_rois
        roi=rois(cnt_roi);
    
        path_out=fullfile(path_home,join(['Model',num2str(model),'_ROI',num2str(roi)]));
        if exist(path_out)~=7
            mkdir(path_out);
        end
        inputs{1, cnt_roi} = {path_out}; % Factorial design specification: Directory - cfg_files
    
        scans={};
        for i=1:length(clinical_data)
            path_zroifc=fullfile(path_home,'FC_FunRawARCWSF',join(['zROI',num2str(roi),'FCMap_CSUB-',num2str(clinical_data(i),'%5.5u'),'C-01.nii,1']));
            scans(i,1)={path_zroifc};
        end
        inputs{2, cnt_roi} = scans; % Factorial design specification: Scans - cfg_files
    end
    
    spm('defaults', 'PET');
    spm_jobman('run', jobs, inputs{:});
    
end
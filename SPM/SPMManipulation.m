% Matlab script to modify spm batch file

path_in = 'H:\MRI\pnTTC\Prosociality_DC_Dr_Okada\SPM\util';
file_in = 'batch_model.mat';
path_out = 'H:\MRI\pnTTC\Prosociality_DC_Dr_Okada\SPM';

%models = [1,2,3];
%rois = [1,2,3,4,5,6,7,8];
models = [1];
rois = [1];
n_covariates=6;


spm_in = load(fullfile(path_in,file_in));
for model=models
    for roi=rois
        spm_out = spm_in;
        %disp(spm_out.matlabbatch{1,1}.spm.stats.factorial_design.dir)
        path_model_roi=fullfile(path_out,join(['Model',num2str(model)]),join(['ROI',num2str(roi)]));
        spm_out.matlabbatch{1,1}.spm.stats.factorial_design.dir={fullfile(path_model_roi,'Stat')};
        %disp(spm_out.matlabbatch{1,1}.spm.stats.factorial_design.dir)
        %id_group=load(fullfile(path_in,join(['id_model',num2str(model),'.txt'])));
        clinical_data=csvread(fullfile(path_in,'model1.csv'),1,0);
        scans={};
        for i=1:length(clinical_data)
            scans(i,1)={fullfile(path_model_roi,'FC',join(['zROI',num2str(roi),'FCMap_CSUB-',num2str(clinical_data(i),'%5.5u'),'C-01.nii,1']))};
        end
        spm_out.matlabbatch{1,1}.spm.stats.factorial_design.des.mreg.scans=scans;
        
        for i=1:n_covariates
            spm_out.matlabbatch{1,1}.spm.stats.factorial_design.des.mreg.mcov(i).c=clinical_data(:,i+1);
        end
        matlabbatch=spm_out.matlabbatch;
        save(fullfile(path_model_roi,'Stat','model.mat'),'matlabbatch');
    end
end
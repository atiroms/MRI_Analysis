% List of open inputs
% Factorial design specification: Directory - cfg_files
% Factorial design specification: Scans - cfg_files
%nrun = X; % enter the number of runs here
jobfile = {'C:\Users\NICT_WS\GitHub\MRI_Analysis\spm\batch\batch3_glm_cocoro\batch03_job.m'};

% original parameters
path_home='D:\MRI_img\dc\dc_re';
list_roi=[37,38,41,42,71,72,73,74,75,76,77,78];
list_group=[1,3];
n_roi=length(list_roi);

jobs = repmat(jobfile, 1, n_roi);
inputs = cell(2, n_roi);
for i_roi = 1:n_roi
    roi=num2str(list_roi(i_roi));
    
    % prepare stat folder (first unspecified parameter)
    path_stat=fullfile(path_home,'stat',join(['roi_',roi]));
    if exist(path_stat)~=7
        mkdir(path_stat);
    end
    
    % prepare input nifti files (second unspecified parameter)
    list_path_img={};
    for group = list_group
        path_dir_img=fullfile(path_home,'zroi',join(['zROI',roi,'_',num2str(group)]));
        
        path_subj=fullfile(path_home,'stat',join(['mrid_',num2str(group),'.csv']));
        file_subj=fopen(path_subj);
        id_subj=fgetl(file_subj);
        while ischar(id_subj)
            file_img=fullfile(path_dir_img,join(['zROI',roi,'FCMap_F_',id_subj,'.nii']));
            list_path_img(length(list_path_img)+1,1)={file_img};
            %disp(id_subj)
            id_subj=fgetl(file_subj);
        end
        fclose(file_subj);
    end
    
    inputs{1, i_roi} = {path_stat}; % Factorial design specification: Directory - cfg_files
    inputs{2, i_roi} = list_path_img; % Factorial design specification: Scans - cfg_files
end
spm('defaults', 'PET');
spm_jobman('run', jobs, inputs{:});

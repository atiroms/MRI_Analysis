% List of open inputs
% Optimal Thresholding: Input Image - cfg_files
nrun = X; % enter the number of runs here
jobfile = {'D:\MRI\pnTTC\c1c2_struc\spm\02_spm\log\batch_preproc_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(1, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Optimal Thresholding: Input Image - cfg_files
end
spm('defaults', 'PET');
spm_jobman('run', jobs, inputs{:});

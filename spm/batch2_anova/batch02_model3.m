% List of open inputs
% Factorial design specification: Directory - cfg_files
% Factorial design specification: Scans - cfg_files
% Factorial design specification: Scans - cfg_files
% Factorial design specification: Scans - cfg_files
% Factorial design specification: Scans - cfg_files
nrun = X; % enter the number of runs here
jobfile = {'/home/atiroms/Documents/GitHub/MRI_Analysis/SPM/batch2_anova/batch02_model3_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(5, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Factorial design specification: Directory - cfg_files
    inputs{2, crun} = MATLAB_CODE_TO_FILL_INPUT; % Factorial design specification: Scans - cfg_files
    inputs{3, crun} = MATLAB_CODE_TO_FILL_INPUT; % Factorial design specification: Scans - cfg_files
    inputs{4, crun} = MATLAB_CODE_TO_FILL_INPUT; % Factorial design specification: Scans - cfg_files
    inputs{5, crun} = MATLAB_CODE_TO_FILL_INPUT; % Factorial design specification: Scans - cfg_files
end
spm('defaults', 'PET');
spm_jobman('run', jobs, inputs{:});

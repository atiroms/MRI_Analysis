% List of open inputs
% Factorial design specification: Scans - cfg_files
nrun = X; % enter the number of runs here
jobfile = {'C:\Users\NICT_WS\GitHub\MRI_Analysis\spm\batch\batch3_glm_cocoro\batch01_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(1, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Factorial design specification: Scans - cfg_files
end
spm('defaults', 'PET');
spm_jobman('run', jobs, inputs{:});

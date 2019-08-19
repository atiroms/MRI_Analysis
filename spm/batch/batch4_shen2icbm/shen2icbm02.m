% List of open inputs
nrun = X; % enter the number of runs here
jobfile = {'C:\Users\NICT_WS\GitHub\MRI_Analysis\spm\batch\batch4_shen2icbm\shen2icbm02_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'PET');
spm_jobman('run', jobs, inputs{:});

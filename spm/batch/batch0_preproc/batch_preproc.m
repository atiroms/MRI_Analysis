% List of open inputs
nrun = X; % enter the number of runs here
jobfile = {'D:\MRI\pnTTC\c1c2_struc\spm\02_spm\log\batch_preproc_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'PET');
spm_jobman('run', jobs, inputs{:});

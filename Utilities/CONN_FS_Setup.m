
% Parameters
TR=2.5;
working_dir='D:\atiroms\MRI\iTTC_rsfMRI_C\CONN_FS\05_preproc';
%working_dir='G:\MRI\iTTC_rsfMRI_C\CONN_FS\05_test4';
%working_dir='G:\MRI\iTTC_rsfMRI_C\CONN_FS\10_test9';
project_filename='conn_project05.mat';
%project_filename='conn_project10.mat';
functional_prefix='\functional\CSUB-';
functional_suffix='C-01.nii';
structural_prefix='\structural\';
structural_suffix='\mri\T1.mgz';
roi_suffix='\mri\aparc+aseg.mgz';


subject_id=importdata('subject_id.txt');
subject_id=subject_id(:,2);
n_subject=size(subject_id,1);

cd(working_dir);

batch.filename=fullfile(pwd,project_filename);

batch.Setup.isnew=1;
batch.Setup.nsubjects=n_subject;
batch.Setup.RT=TR;

for i=1:n_subject
    idstr=num2str(subject_id(i,1),'%5.5u');
    functional_filename=strcat(working_dir,functional_prefix,idstr,functional_suffix);
    structural_filename=strcat(working_dir,structural_prefix,idstr,structural_suffix);
    batch.Setup.functionals{i}{1}=functional_filename;
    batch.Setup.structurals{i}=structural_filename;
end

batch.Setup.rois.names={'fs'};
%batch.Setup.rois.names={'Grey Matter','White Matter','CSF','fs'};

for i=1:n_subject
    idstr=num2str(subject_id(i,1),'%5.5u');
    roi_filename=strcat(working_dir,structural_prefix,idstr, roi_suffix);
    batch.Setup.rois.files{1}{i}=roi_filename;
end

conn_batch(batch);

conn
function conn_batch_MR3(varargin)
%
% batch processing script for rs-fMRI
% 
% This script is based on conn_batch_NYU.m
%
% Usage:
% Run conn_batch_MR3.m
% The script will:
%    a) sort functional and structure files into subdirectories.
%    b) Preprocessing of the anatomical and functional volumes
%       (normalization & segmentation of anatomical volumes; realignment,
%        coregistration, normalization, outlier detection, 
%	     and smooting of the functional volumes)
%    c) Estimate first-level seed-to-voxel connectivity maps for each of 
%       the default seeds (located in the conn/rois folder), separately 
%       for each subject and for each of the three test-retest sessions.
% Kiyotaka Nemoto 15 Oct 2016

% Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these variable first %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TR=2.5; % Repetition time = 2.5 seconds
sorder='ascending';
    %Select one from below
    % 'ascending'
    % 'descending'
    % 'interleaved (middle-top)'
    % 'interleaved (bottom-up)'
    % 'interleaved (top-down)'
    % 'interleaved (Siemens)'  
nremove=10;          %initial scans to remove
t1vprefix='V_';
rsfprefix='F_';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Display settings
display(['TR: ' num2str(TR) 's']);
display(['Slice Order: ' sorder]);
display(['Initial scans to be romoved: ' num2str(nremove)]);

%% Sort files into subdirectories
imglist=spm_select(Inf,'image','Specify functional images',{},pwd,rsfprefix);
[cwd]=fileparts(imglist(1,:));
cd(cwd); %change directory to the one in which image files exist.
display('Sort files into subdirectories');

for i=1:size(imglist,1)
    func_image=imglist(i,:);
    [path func_image_file ext]=fileparts(func_image);
    ID=func_image_file(:,3:end);
    dircheck=exist(ID,'dir');
    if dircheck~=7
        mkdir(ID);
    end
    t1volume=[t1vprefix ID '.nii'];
    rsfmri=[rsfprefix ID '.nii'];
    dircheck=exist('original','dir');
    if dircheck~=7
        mkdir('original');
    end
    dircheck=exist('results','dir');
    if dircheck~=7
        mkdir('results')
    end
    copyfile(t1volume,'original');
    copyfile(rsfmri,'original');
    movefile(t1volume,ID);
    movefile(rsfmri,ID);
    
    %% Identify functional/structural files
    cd(ID);
    FUNCTIONAL_FILE=cellstr(conn_dir(rsfmri));
    STRUCTURAL_FILE=cellstr(conn_dir(t1volume));

    %% CONN-SPECIFIC SECTION: 
    %% RUNS PREPROCESSING/SETUP/DENOISING/ANALYSIS STEPS
    %% Prepares batch structure
    clear batch;
    projectfile=['conn_' ID '.mat'];
    batch.filename=fullfile(cwd,ID,projectfile); % New project name
    
    %% SETUP & PREPROCESSING step
    %% (using default values for most parameters, see help conn_batch 
    %% to define non-default values)
    % CONN Setup
	% Default options
	% (uses all ROIs in conn/rois/ directory); 
	% see conn_batch for additional options 
    % CONN Setup.preprocessing		
	% (realignment/coregistration/segmentation/normalization/smoothing)
    batch.Setup.isnew=1;
    batch.Setup.nsubjects=1;
    batch.Setup.RT=TR;
		% TR (seconds)
    batch.Setup.functionals=repmat({{}},[1,1]);
		% Point to functional volumes for each subject/session
    batch.Setup.functionals{1}{1}{1}=FUNCTIONAL_FILE{1,1}; 
    batch.Setup.structurals=STRUCTURAL_FILE;
		% Point to anatomical volumes for each subject
    %batch.Setup.rois.files=
    batch.Setup.conditions.names={'rest'};
    batch.Setup.conditions.durations{1}{1}{1}=inf;
    batch.Setup.preprocessing.steps={'functional_removescans',...
        'functional_realign&unwarp','functional_center',...
	'functional_slicetime','structural_center',...
	'structural_segment&normalize','functional_normalize',...
        'functional_art','functional_smooth'};
    batch.Setup.preprocessing.removescans=nremove;
    batch.Setup.preprocessing.sliceorder=sorder;
    batch.Setup.preprocessing.art_thresholds(1)=3; %conservative
    batch.Setup.preprocessing.art_thresholds(2)=0.5; %conservative
    batch.Setup.done=1;
    batch.Setup.overwrite='Yes';                            
        
    %% DENOISING step
    % CONN Denoising
	% Default options 
	% (uses White Matter+CSF+realignment+scrubbing+
	%  conditions as confound regressors); 
	% see conn_batch for additional options 
    batch.Denoising.filter=[0.01, 0.1];
		% frequency filter (band-pass values, in Hz)
    batch.Denoising.done=1;
    batch.Denoising.overwrite='Yes';
    
    %% FIRST-LEVEL ANALYSIS step
    % CONN Analysis
	% Default options
	% (uses all ROIs in conn/rois/ as connectivity sources); 
	% see conn_batch for additional options 
    batch.Analysis.done=1;
    batch.Analysis.overwrite='Yes';
    
    %% Run all analyses
    conn_batch(batch);
    
    %% CONN Display
    % launches conn gui to explore results
    %conn
    %conn('load',fullfile(cwd,projectfile));
    %conn gui_results
    
    %% Extract correlation matrix in DMN
    projectdir=['conn_' ID];
    resultmat=fullfile(cwd,ID,projectdir,'results','firstlevel','ANALYSIS_01','resultsROI_Subject001_Condition001.mat');
    load(resultmat);
    column_name=names(1,:);
    corrmatrix=Z(:,1:end-1);

    %% generate csv
    resultfilename=[ID '_corr_matrix.csv'];
    resultfile=fullfile(cwd,'results',resultfilename);
    csvwrite(resultfile,corrmatrix);
    
    roinamefile=fullfile(cwd,'results','roinames.csv');
    fid=fopen(roinamefile,'wt');
    fprintf(fid,'%s\n',column_name{:});
    fclose(fid);
    
    %% Save figure of correlation matrix
    figure;
    colormap(jet);
    imagesc(corrmatrix);
    axis square;
    caxis([-0.5 1.0]);
    colorbar;
    pngfilename=[ID '_corr_matrix.png'];
    pngfile=fullfile(cwd,'results',pngfilename);
    saveas(gcf,pngfile);
    close(gcf);
    
    %% Save result of scrubbing
    artfilename=[ID '_scrubbing.jpg'];
    copyfile('art_screenshot.jpg',artfilename);
    movefile(artfilename,'../results');
    
    %% change to the upper directory
    cd ../;
end

display('Done');


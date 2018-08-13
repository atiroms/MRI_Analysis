% specify radius of spheres to build in mm
radiusmm = 5;

coordinates=load('PowerAtlas.txt');
% Specify Output Folders for two sets of images (.img format and .mat format)
roi_dir_img = 'img';
%roi_dir_mat = 'mat';
% Make an img and an mat directory to save resulting ROIs
mkdir('img');
%mkdir('mat');
% Go through each set of coordinates from the specified file (line 2)
spherelistrows = length(coordinates(:,1));
for spherenumbers = 1:spherelistrows
    % maximum is specified as the centre of the sphere in mm in MNI space
    maximum = coordinates(spherenumbers,1:3);
    sphere_centre = maximum;
    sphere_radius = radiusmm;
    sphere_roi = maroi_sphere(struct('centre', sphere_centre,'radius', sphere_radius));

    % Define sphere name using coordinates
    coordsx = num2str(maximum(1));
    coordsy = num2str(maximum(2));
    coordsz = num2str(maximum(3));
    spherelabel = sprintf('%s_%s_%s', coordsx, coordsy, coordsz);
    sphere_roi = label(sphere_roi, spherelabel);

    % save ROI as MarsBaR ROI file

    %matfilename=sprintf('%3.3u_%dmm_%s_roi.mat', spherenumbers, radiusmm, spherelabel);
    %saveroi(sphere_roi, fullfile(roi_dir_mat, matfilename));
    % Save as image
    sp = maroi('classdata', 'spacebase');
    niifilename=sprintf('%3.3u_%dmm_%s_roi.nii', spherenumbers, radiusmm, spherelabel);
    save_as_image(sphere_roi, fullfile(roi_dir_img, niifilename),sp);

end
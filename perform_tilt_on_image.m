function [ im_out ] = perform_tilt_on_image( im1,t, theta )
%PERFORM_IMAGE_TILT - This function compiles if needed the Mex-Function
%                     perform_tilt.cpp, and if it's already compiled just
%                     calls that function!

% NOTES:
% Matlab handles differently images with respect to the c++ code!
% We need to use the image transpose instead

ALWAYS_COMPILE = false;

currentfolder = pwd; % this saves the current folder
cd ./simu_tilts_cpp/;

if ( (exist('perform_tilt')==0)||ALWAYS_COMPILE )
    mex -v -I. perform_tilt.cpp frot.cpp splines.cpp fproj.cpp library.cpp ...
        flimage.cpp filter.cpp numerics1.cpp
    pause(0.2);    
end

% Calling the Mex-Function
if ( exist('perform_tilt')==3)
    im_out = perform_tilt(im1',t,theta);
    im_out = im_out';
else
    im_out = 0;
    error('Error: The Mex-Function perform_tilt was not compiled.')
end

cd(currentfolder); % go back to the main directory
end


function [ x_out,y_out, inside ] = get_tilted_coor( x, y, w, h, t, theta, direction)
%PERFORM_IMAGE_TILT - This function compiles if needed the Mex-Function
%                     perform_tilt.cpp, and if it's already compiled just
%                     calls that function!

% NOTES:
% Matlab handles differently images with respect to the c++ code!
% We need to use the image transpose instead

ALWAYS_COMPILE = false;
BORDER_SAFETY = 0; % only for the inverse transform

currentfolder = pwd; % this saves the current folder
cd ./simu_tilts_cpp/;

if ( (exist('perform_tilt')==0)||ALWAYS_COMPILE )
    mex -v -I. tilted_coor.cpp numerics1.cpp libNumerics/numerics.cpp
    pause(0.2);    
end

% Calling the Mex-Function
if ( exist('tilted_coor')==3)
    data_out = tilted_coor(x, y, w, h, t, theta, direction,BORDER_SAFETY);
    x_out = data_out(1);
    y_out = data_out(2);
    inside = (data_out(3)>0);
else
    error('Error: The Mex-Function perform_tilt was not compiled.')
end

cd(currentfolder); % go back to the main directory
end


function [ data_stats ] = get_sift_stats( im1)
%PERFORM_SIFT Summary of this function goes here

% NOTES:
% Matlab handles differently images with respect to the c++ code!
% We need to use the image transpose instead

ALWAYS_COMPILE = true;


setenv('OMP_NUM_THREADS', '4'); % use 4 threads for computing
currentfolder = pwd; % this saves the current folder
cd ./sift/;
if ( (exist('get_sift_stats_mex')==0) || ALWAYS_COMPILE )
        mex -v -I. -I"/home/rdguez-mariano/Sources/opencv_3.2.0/include/opencv" ...
            -I"/home/rdguez-mariano/Sources/opencv_3.2.0/include" ...
            -I"/home/rdguez-mariano/Sources/opencv_3.2.0/share/OpenCV" ...
            -I"/home/rdguez-mariano/Sources/opencv_3.2.0/include/opencv2" ...
            -L"/home/rdguez-mariano/Sources/opencv_3.2.0/lib" ...            
            get_sift_stats_mex.cpp ...
            CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"...
            CFLAGS="\$CFLAGS -fopenmp"...
            -lopencv_legacy -lopencv_imgproc -lopencv_core -lopencv_contrib -lopencv_ml ...
            -lopencv_objdetect -lopencv_calib3d -lopencv_flann -lopencv_features2d -lopencv_video ...
            -lopencv_gpu -lopencv_xfeatures2d;
        
    pause(0.2);
end

% Calling the Mex-Function
if ( exist('get_sift_stats_mex')==3)
    data_stats = get_sift_stats_mex(double(im1'));
    cd(currentfolder); % go back to the main directory
else
    error('Error: The Mex-Function ASIFT_matlab was not compiled.')
end

end

project(get_sift_stats)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)



FIND_PACKAGE(OpenMP)
if (OPENMP_FOUND)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif (OPENMP_FOUND)


# Not the MATLAB compiler!!!
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D _NOMEX=1")

set (OpenCV_DIR "/home/rdguez-mariano/Sources/opencv_3.2.0/share/OpenCV")
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/home/rdguez-mariano/Sources/opencv_3.2.0/share/OpenCV")
find_package( OpenCV 3.2 REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )

####### Different Configs For IMAS
set(src
get_sift_stats_mex.cpp
)


add_executable(get_sift_stats_mex ${src})
TARGET_LINK_LIBRARIES(get_sift_stats_mex ${OpenCV_LIBS})

message("***************************************")
message("Include Dirs: \n ${OpenCV_INCLUDE_DIRS}")
message(" ")
message("Libs: \n ${OpenCV_LIBS}")
message("***************************************")

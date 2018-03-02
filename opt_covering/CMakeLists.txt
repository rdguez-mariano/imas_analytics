project(opt_covering)
include(FindOpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    set(LINKER_FLAGS "${LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(gprof OFF)
if (gprof)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
     message("************* intended for gprof *************")
endif()


set(saveflag OFF)
if (saveflag)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D _saveflag=1")
    message("************* SAVEFLAG ON *************")
endif()


cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)


set(PNG ON)
if (PNG)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D _PNG=1")
    message("************* PNG ON *************")
    add_executable(${PROJECT_NAME} ${SRC_LIST} io_png/io_png.c)

    include_directories(.
                                            ./io_png
                                            ./io_png/libs/png)

    TARGET_LINK_LIBRARIES(${PROJECT_NAME} png)

else()
    add_executable(${PROJECT_NAME} ${SRC_LIST})

    include_directories(.)

    TARGET_LINK_LIBRARIES(${PROJECT_NAME})
endif()

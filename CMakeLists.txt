cmake_minimum_required(VERSION 3.7)
project(SiftySifty)

FIND_PACKAGE(OpenCV REQUIRED)

# for speed test should add openmp support
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -lstdc++ -O3")
# set(CMAKE_CXX_COMPILER /usr/local/bin/gcc-7)

FIND_PACKAGE(OpenCV REQUIRED)

include_directories(include)

set(SOURCE_FILES
        main.cpp
        include/utils.h
        include/structs.h
        src/linearfilter.cpp
        include/linearfilter.h
        include/boxfilter.h
        src/boxfilter.cpp
        src/gaussfiler.cpp
        include/gaussfiler.h
        include/siftysifty.h
        src/siftysifty.cpp
        src/iirfilter.cpp
        include/iirfilter.h
        include/filter.h
        src/siftysiftytest.cpp
        include/siftysiftytest.h
        include/imageutils.h)

add_executable(SiftySifty ${SOURCE_FILES})

target_link_libraries(SiftySifty ${OpenCV_LIBS})

# CMakeList.txt : CMake project for SCAN2CAD, include source and define
# project specific logic here.
#
cmake_minimum_required(VERSION 3.8)
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project ("SCAN2CAD")

find_package(PCL 1.2 REQUIRED)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})


add_executable (SCAN2CAD "SCAN2CAD.cpp" "SCAN2CAD.h")
target_link_libraries (SCAN2CAD ${PCL_LIBRARIES})

#add_definitions(-Wall -O3 -DNDEBUG -pipe -ffast-math -funroll-loops -ftree-vectorize -fomit-frame-pointer)
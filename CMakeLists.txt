project(particles10)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
#add_executable(${PROJECT_NAME} ${SRC_LIST})
add_executable(${PROJECT_NAME} main.cpp psystem.cpp)
add_executable(psdistance psdistance.cpp psystem.cpp psystem.h ./intersections/intersections.cpp ./intersections/intersections.h)
set(CMAKE_CXX_FLAGS "-O3 -std=c++11 -Wall -fexceptions -fexpensive-optimizations")


cmake_minimum_required(VERSION 3.1)

project(Parallel_Lab_1)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O0")


add_executable(main main.cpp)
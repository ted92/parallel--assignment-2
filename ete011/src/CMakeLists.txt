cmake_minimum_required(VERSION 3.3)
project(src)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    Makefile
    tsp-par.c
    tsp-seq.c)

add_executable(src ${SOURCE_FILES})
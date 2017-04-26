cmake_minimum_required(VERSION 3.7)
project(projetMatriceAnthonyFusco)

set(CMAKE_C_STANDARD 99)

set(SOURCE_FILES main.c)

include_directories(/usr/include/mpi/)

SET(CMAKE_C_COMPILER mpicc)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp -Wall")

add_executable(projetMatriceAnthonyFusco ${SOURCE_FILES})
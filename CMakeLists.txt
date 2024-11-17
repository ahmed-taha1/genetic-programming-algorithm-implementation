cmake_minimum_required(VERSION 3.2)
project(ass1)

set(CMAKE_CXX_STANDARD 17)

add_executable(ass1 main.cpp
        GeneticAlgorithm.h
        Chromosom.h
        GeneticAlgorithm.cpp
        RandomNumGeneratorHelper.h)

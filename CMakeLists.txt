cmake_minimum_required(VERSION 3.2.0)

project(my-simulation)

find_package(BioDynaMo REQUIRED)
include(${BDM_USE_FILE})
include_directories("src")

bdm_add_executable(my-simulation
                   SOURCES src/my-simulation.cc
                   HEADERS src/my-simulation.h
                   LIBRARIES ${BDM_REQUIRED_LIBRARIES})

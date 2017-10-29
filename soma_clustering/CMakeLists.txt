cmake_minimum_required(VERSION 3.2.0)

project(soma-clustering)

find_package(BioDynaMo REQUIRED)
include(${BDM_USE_FILE})
include_directories("src")

bdm_add_executable(soma_clustering
                   SOURCES src/soma_clustering.cc
                   HEADERS src/soma_clustering.h
                   LIBRARIES ${BDM_REQUIRED_LIBRARIES})

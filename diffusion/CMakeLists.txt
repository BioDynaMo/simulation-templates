cmake_minimum_required(VERSION 3.2.0)

project(diffusion)

find_package(BioDynaMo REQUIRED)
include(${BDM_USE_FILE})
include_directories("src")

bdm_add_executable(diffusion
                   SOURCES src/diffusion.cc
                   HEADERS src/diffusion.h src/diffusion_biology_modules.h
                   LIBRARIES ${BDM_REQUIRED_LIBRARIES})

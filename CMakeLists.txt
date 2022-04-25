cmake_minimum_required(VERSION 3.10)

project(filter VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(EIGEN 3.3 REQUIRED)

add_subdirectory(src)

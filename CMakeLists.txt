cmake_minimum_required(VERSION 3.10)

project(filter VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(EIGEN 3.3 REQUIRED)

# Link in GTest
find_package(GTEST 1.10 REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_subdirectory(src)

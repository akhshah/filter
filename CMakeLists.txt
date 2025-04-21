cmake_minimum_required(VERSION 3.15...3.31)

set(CMAKE_C_COMPILER "gcc-14")
set(CMAKE_CXX_COMPILER "g++-14")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(Filter VERSION 0.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/share/"
    "/opt/homebrew/include")
find_package(Eigen3 REQUIRED NO_MODULE)

add_executable("${PROJECT_NAME}" ./main.cpp)

target_link_libraries("${PROJECT_NAME}" Eigen3::Eigen)

cmake_minimum_required(VERSION 3.10)
project(finiteVolumeSolver)

set(CMAKE_CXX_STANDARD 14)
add_subdirectory(Grids)
find_package(Eigen3 3.3)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin")
add_executable(finiteVolumeSolver main.cpp Grids/Grid.cpp Grids/Grid.h Grids/CartesianGrid.cpp Grids/CartesianGrid.h)
target_link_libraries(finiteVolumeSolver Eigen3::Eigen)
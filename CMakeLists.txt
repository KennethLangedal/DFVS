cmake_minimum_required(VERSION 3.0.0)
project(dfvs VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED)
# set(CMAKE_CXX_FLAGS "-D VERBOSE")
set(CMAKE_CXX_FLAGS "-O3")

include_directories(include/)

# add_executable(dfvs main.cpp src/bitvector.cpp src/graph.cpp src/reductions.cpp src/solver.cpp src/local_search.cpp)

add_executable(dfvs_heuristic main_heuristic.cpp src/sparse_graph.cpp src/sparse_local_search.cpp src/bitvector.cpp src/sparse_reductions.cpp)
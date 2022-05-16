cmake_minimum_required(VERSION 3.0.0)
project(dfvs VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED)
set(CMAKE_CXX_FLAGS "-D VERBOSE")
# set(CMAKE_CXX_FLAGS "-O3 -D NDEBUG")

include_directories(include/)

add_executable(
  dfvs_exact
  main_exact.cpp
  src/sparse_graph.cpp
  src/bitvector.cpp
  src/reductions.cpp
  src/reduction_engine.cpp
  src/local_search.cpp
  src/solver.cpp
  src/bounds.cpp)

add_executable(
  dfvs_heuristic
  main_heuristic.cpp
  src/sparse_graph.cpp
  src/bitvector.cpp
  src/reductions.cpp
  src/reduction_engine.cpp
  src/local_search.cpp
  src/local_search_edgew.cpp
  src/bounds.cpp)

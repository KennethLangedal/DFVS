cmake_minimum_required(VERSION 3.0.0)
project(dfvs VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED)

include_directories(include/)

include(CTest)
enable_testing()

add_executable(dfvs main.cpp src/graph.cpp src/reductions.cpp src/solver.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

# file(GLOB files "data/exact/*")
# foreach(file ${files})
#     message(STATUS "X=${file}")
# endforeach()
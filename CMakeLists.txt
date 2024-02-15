cmake_minimum_required(VERSION 3.25)
project(ModularDecomposition_v3)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -DNDEBUG -march=native -mtune=native")

include_directories(ModularDecomposition_v3)

file(GLOB sources ModularDecomposition_v3/*.cpp)

add_executable(ModularDecomposition ${sources})

cmake_minimum_required(VERSION 2.6)

project(smg2s-test)

add_executable(smg2s-test test.cpp)
target_link_libraries(smg2s-test parvec-lib parmat-lib)


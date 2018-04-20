#!/bin/bash

rm -rf CMakeCache.txt CMakeFiles cmake_install.cmake 

find . -maxdepth 3 -type d \( ! -name . \) -exec bash -c "cd '{}' && rm -rf CMakeCache.txt CMakeFiles cmake_install.cmake" \;



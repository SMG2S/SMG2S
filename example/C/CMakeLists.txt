add_executable(c_interface.exe test_c.c)
target_link_libraries(c_interface.exe PRIVATE MPI::MPI_CXX smg2s2c)
target_include_directories(c_interface.exe PRIVATE ${CMAKE_SOURCE_DIR})


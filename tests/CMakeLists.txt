aux_source_directory(. DIR_SRC)

add_executable(test_nilpotent.exe test_nilpotent.cpp)


add_executable(test_parVector.exe test_parVector.cpp)
target_link_libraries(test_parVector.exe PRIVATE PRIVATE MPI::MPI_CXX)
target_include_directories(test_parVector.exe PRIVATE ${MPI_CXX_INCLUDE_PATH})


add_executable(test_parMatrix.exe test_parMatrix.cpp)
target_link_libraries(test_parMatrix.exe PRIVATE PRIVATE MPI::MPI_CXX)
target_include_directories(test_parMatrix.exe PRIVATE ${MPI_CXX_INCLUDE_PATH})
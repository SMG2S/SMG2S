aux_source_directory(. DIR_SRC)

add_executable(run.exe test.cpp)

link_directories(${SMG2S_SOURCE_FILE}/${DIR_PARVEC_SRC} ${SMG2S_SOURCE_FILE}/${DIR_PARMAT_SRC})

target_link_libraries(run.exe PRIVATE ${MPI_C_LIBRARIES} parvec parmat)

target_include_directories(run.exe PRIVATE ${MPI_C_INCLUDE_PATH})


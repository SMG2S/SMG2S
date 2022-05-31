add_executable(ex1.exe ex1.cpp)
target_link_libraries(ex1.exe PRIVATE smg2s)

add_executable(ex2.exe ex2.c)
target_link_libraries(ex2.exe PRIVATE smg2s2c)

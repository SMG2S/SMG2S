add_executable(test_nilpotent.exe test_nilpotent.cpp)
target_link_libraries(test_nilpotent.exe PRIVATE smg2s )

add_executable(test_parVector.exe test_parVector.cpp)
target_link_libraries(test_parVector.exe PRIVATE smg2s)

add_executable(test_parMatrix.exe test_parMatrix.cpp)
target_link_libraries(test_parMatrix.exe PRIVATE smg2s)

add_executable(test_smg2s.exe test_smg2s.cpp)
target_link_libraries(test_smg2s.exe PRIVATE smg2s)

add_executable(test_smg2s_c.exe test_smg2s_c.c)
target_link_libraries(test_smg2s_c.exe PRIVATE smg2s2c)
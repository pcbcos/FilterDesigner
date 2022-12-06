cmake_minimum_required(VERSION 3.23)
project(EllipticFilter CXX C)


set(CMAKE_CXX_STANDARD 23)

include_directories("P:/LIBS/GMP_fast/include")
include_directories("P:/LIBS/mpfr/include")
include_directories(include)
aux_source_directory(src src)

link_directories(P:/LIBS/GMP_fast/bin)
link_directories(P:/LIBS/mpfr/bin)

add_executable(AFD_test AFD_test.cpp ${src})
add_executable(DFD_test DFD_test.cpp ${src})
#add_executable(test main2.cpp)
#add_executable(test2 main3.c)
target_link_libraries(AFD_test libgmp-10.dll libgmpxx-4.dll libmpfr-6.dll)
target_link_libraries(DFD_test libgmp-10.dll libgmpxx-4.dll libmpfr-6.dll)
#target_link_libraries(test libgmp-10.dll libgmpxx-4.dll libmpfr-6.dll)
cmake_minimum_required(VERSION 3.23)
project(EllipticFilter CXX C)


set(CMAKE_CXX_STANDARD 23)

include_directories("thirdpart/GMP/include")
include_directories("thirdpart/mpfr/include")
include_directories(include)
aux_source_directory(src src)

link_directories("thirdpart/GMP/lib")
link_directories("thirdpart/mpfr/lib")

add_executable(AFD_test AFD_test.cpp ${src})
add_executable(DFD_test DFD_test.cpp ${src})
#add_executable(test main2.cpp)
#add_executable(test test.cpp)
target_link_libraries(AFD_test  mpfr.a gmp.a gmpxx.a)
target_link_libraries(DFD_test libgmp.a libgmpxx.a libmpfr.a)
#target_link_libraries(test libgmp-10.dll libgmpxx-4.dll libmpfr-6.dll)
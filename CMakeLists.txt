cmake_minimum_required(VERSION 2.8.12)
project("TI Uniformity")
#set(CMAKE_CXX_FLAGS "-std=c++11 -g")
set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -pg -fprofile-arcs -ftest-coverage")
set(PROGRAM_NAME "uniformity.out")

set(SRC "main.cpp")

add_executable(${PROGRAM_NAME} ${SRC})
target_link_libraries(${PROGRAM_NAME} pthread)

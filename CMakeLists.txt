cmake_minimum_required(VERSION 3.29)
project(AIS4104-assignment-n)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(src)
#add_subdirectory(screws)
add_subdirectory(math)
add_dependencies(mainsrc math)

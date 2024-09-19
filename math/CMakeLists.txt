find_package(Eigen3 CONFIG REQUIRED)

add_library(math math.cpp
        include/math/math.h)

target_link_libraries(math PRIVATE Eigen3::Eigen)

target_include_directories(math PUBLIC include/)
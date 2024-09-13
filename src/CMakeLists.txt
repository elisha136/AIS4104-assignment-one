find_package(Eigen3 CONFIG REQUIRED)

add_executable(mainsrc main.cpp)

target_link_libraries(mainsrc PUBLIC math PRIVATE Eigen3::Eigen)

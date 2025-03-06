cmake_minimum_required(VERSION 3.10)
project(Cubic3D)

add_executable(cubic3d cubic3d.cpp)  # Создать исполняемый файл

target_link_libraries(cubic3d ${VTK_LIBRARIES})  # Прилинковать библиотеки VTK
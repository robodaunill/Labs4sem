cmake_minimum_required(VERSION 3.10)
project(Cubic2D)

add_executable(cubic2d cubic2d.cpp)  # Создать исполняемый файл

target_link_libraries(cubic2d ${VTK_LIBRARIES})  # Прилинковать библиотеки VTK
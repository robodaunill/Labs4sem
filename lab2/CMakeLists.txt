cmake_minimum_required(VERSION 3.10)  # Или выше, в зависимости от требований VTK
project(VTKExamples)

find_package(VTK REQUIRED)  # Найти установленный VTK
include(${VTK_USE_FILE})    # Подключить настройки VTK

add_subdirectory(cubic2d)
add_subdirectory(cubic3d)
add_subdirectory(tetr3d)
#include <gmsh.h>

int main(int argc, char **argv)
{
    // Инициализация Gmsh
    gmsh::initialize();

    // Создание новой модели
    gmsh::model::add("torus");

    // Параметры тора
    double R = 5.0;  // Большой радиус тора
    double r = 1.0;  // Малый радиус тора (толщина стенки)

    // Создание тора
    int tag = gmsh::model::occ::addTorus(0, 0, 0, R, r);
    gmsh::model::occ::synchronize();

    // Задание размеров сетки
    double meshSize = r / 4.0; // Размер элемента сетки (чтобы было 3-4 тетраэдра на толщину стенки)
    gmsh::option::setNumber("Mesh.CharacteristicLengthMin", meshSize);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", meshSize);

    // Генерация сетки
    gmsh::model::mesh::generate(3); // 3D тетраэдральная сетка

    // Сохранение сетки в файл
    gmsh::write("torus_mesh.msh");

    // Запуск GUI Gmsh для визуализации
    gmsh::fltk::run();

    // Завершение работы Gmsh
    gmsh::finalize();

    return 0;
}
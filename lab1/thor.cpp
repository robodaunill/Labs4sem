#include <gmsh.h>

int main(int argc, char **argv)
{
    // Инициализация Gmsh
    gmsh::initialize();

    // Создание новой модели
    gmsh::model::add("stl_model");

    // Импорт STL-файла
    gmsh::merge("/home/robodaniil/Labs4sem/lab1/cat-v3.stl"); 

    // Создание поверхности из импортированного STL
    gmsh::model::occ::synchronize();

    // Получаем все поверхности (2D-объекты)
    std::vector<std::pair<int, int>> surfaces;
    gmsh::model::getEntities(surfaces, 2);

    // Создаём объём из поверхностей
    if (surfaces.size() > 0) {
        // Собираем теги поверхностей в вектор
        std::vector<int> surfaceTags;
        for (const auto &surface : surfaces) {
            surfaceTags.push_back(surface.second); // Используем второй элемент пары (тег поверхности)
        }

        // Создаём объём
        int volumeTag;
        gmsh::model::occ::addVolume(surfaceTags, volumeTag);
        gmsh::model::occ::synchronize();
    }

    // Задание размеров сетки
    double meshSize = 1.0; // Размер элемента сетки (можно настроить под ваши нужды)
    gmsh::option::setNumber("Mesh.CharacteristicLengthMin", meshSize);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", meshSize);

    // Генерация сетки
    gmsh::model::mesh::generate(3); // 3D тетраэдральная сетка

    // Сохранение сетки в файл
    gmsh::write("stl_mesh.msh");

    // Запуск GUI Gmsh для визуализации
    gmsh::fltk::run();

    // Завершение работы Gmsh
    gmsh::finalize();

    return 0;
}
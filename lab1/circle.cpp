#include <iostream>
#include <set>
#include <gmsh.h>
#include <cmath>

int main(int argc, char **argv) {
    gmsh::initialize();
    gmsh::model::add("circle");

    double lc = 0.1;
    double radius = 1.0;
    int num_segments = 36;  // Количество сегментов для аппроксимации круга

    std::vector<int> points;
    std::vector<int> arcs;

    // Создаем точку в центре окружности
    int center_tag = gmsh::model::geo::addPoint(0, 0, 0, lc);

    // Создаем точки на окружности
    for (int i = 0; i < num_segments; ++i) {
        double angle = 2 * M_PI * i / num_segments;
        double x = radius * cos(angle);
        double y = radius * sin(angle);
        int point_tag = gmsh::model::geo::addPoint(x, y, 0, lc);
        points.push_back(point_tag);
    }

    // Создаем дуги, соединяющие точки
    for (size_t i = 0; i < points.size(); ++i) {
        int start_point = points[i];
        int end_point = points[(i + 1) % points.size()];
        int arc_tag = gmsh::model::geo::addCircleArc(start_point, center_tag, end_point);
        arcs.push_back(arc_tag);
    }

    // Создаем контур (замкнутую кривую) из дуг
    int curve_loop_tag = 1;
    gmsh::model::geo::addCurveLoop(arcs, curve_loop_tag);

    // Создаем поверхность (круг), ограниченную контуром
    int surface_tag = 1;
    gmsh::model::geo::addPlaneSurface({curve_loop_tag}, surface_tag);

    gmsh::model::geo::synchronize();

    // Добавляем физическую группу (обязательно для сохранения)
    std::vector<int> all_surfaces = {1};
    gmsh::model::addPhysicalGroup(2, all_surfaces, 1);
    gmsh::model::setPhysicalName(2, 1, "Circle");

    gmsh::model::mesh::generate(2);
    gmsh::write("circle.msh");

    std::set<std::string> args(argv, argv + argc);
    if (!args.count("-nopopup")) gmsh::fltk::run();

    gmsh::finalize();
    return 0;
}
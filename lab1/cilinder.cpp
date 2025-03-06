#include <iostream>
#include <set>
#include <gmsh.h>
#include <cmath>

int main(int argc, char **argv) {
    gmsh::initialize();
    gmsh::model::add("cylinder");

    double lc = 0.1;
    double radius = 1.0;
    double height = 2.0;
    int num_segments = 36;  // Количество сегментов для аппроксимации круга

    std::vector<int> bottom_points;
    std::vector<int> bottom_arcs;
    std::vector<int> top_points;
    std::vector<int> top_arcs;

    // Создаем точки в центре нижнего и верхнего кругов
    int bottom_center_tag = gmsh::model::geo::addPoint(0, 0, 0, lc);
    int top_center_tag = gmsh::model::geo::addPoint(0, 0, height, lc);

    // Создаем точки для нижнего круга
    for (int i = 0; i < num_segments; ++i) {
        double angle = 2 * M_PI * i / num_segments;
        double x = radius * cos(angle);
        double y = radius * sin(angle);
        int point_tag = gmsh::model::geo::addPoint(x, y, 0, lc);
        bottom_points.push_back(point_tag);
    }

    // Создаем дуги для нижнего круга
    for (size_t i = 0; i < bottom_points.size(); ++i) {
        int start_point = bottom_points[i];
        int end_point = bottom_points[(i + 1) % bottom_points.size()];
        int arc_tag = gmsh::model::geo::addCircleArc(start_point, bottom_center_tag, end_point);
        bottom_arcs.push_back(arc_tag);
    }

    // Создаем точки для верхнего круга
    for (int i = 0; i < num_segments; ++i) {
        double angle = 2 * M_PI * i / num_segments;
        double x = radius * cos(angle);
        double y = radius * sin(angle);
        int point_tag = gmsh::model::geo::addPoint(x, y, height, lc);
        top_points.push_back(point_tag);
    }

    // Создаем дуги для верхнего круга
    for (size_t i = 0; i < top_points.size(); ++i) {
        int start_point = top_points[i];
        int end_point = top_points[(i + 1) % top_points.size()];
        int arc_tag = gmsh::model::geo::addCircleArc(start_point, top_center_tag, end_point);
        top_arcs.push_back(arc_tag);
    }

    // Создаем контуры (замкнутые кривые) из дуг
    int bottom_curve_loop_tag = 1;
    gmsh::model::geo::addCurveLoop(bottom_arcs, bottom_curve_loop_tag);

    int top_curve_loop_tag = 2;
    gmsh::model::geo::addCurveLoop(top_arcs, top_curve_loop_tag);

    // Создаем поверхности (круги), ограничивающие цилиндр
    int bottom_surface_tag = 1;
    gmsh::model::geo::addPlaneSurface({bottom_curve_loop_tag}, bottom_surface_tag);

    int top_surface_tag = 2;
    gmsh::model::geo::addPlaneSurface({top_curve_loop_tag}, top_surface_tag);

    // Создаем боковые линии, соединяющие круги
    std::vector<int> side_lines;
    for (int i = 0; i < num_segments; ++i) {
        int line_tag = gmsh::model::geo::addLine(bottom_points[i], top_points[i]);
        side_lines.push_back(line_tag);
    }

    //Боковая поверхность цилиндра
    std::vector<int> side_curve_loop_tags;
    for(int i = 0; i < num_segments; ++i) {
        int curve_loop_tag = 100 + i; //уникальный тег
        std::vector<int> curve_loop = {bottom_arcs[i], side_lines[(i+1) % num_segments], -top_arcs[i], -side_lines[i]};
        gmsh::model::geo::addCurveLoop(curve_loop, curve_loop_tag);
        side_curve_loop_tags.push_back(curve_loop_tag);
    }

     std::vector<int> side_surface_tags;
     for(int i = 0; i < num_segments; ++i){
        int surface_tag = 100 + i;
        gmsh::model::geo::addPlaneSurface({side_curve_loop_tags[i]}, surface_tag);
        side_surface_tags.push_back(surface_tag);
     }


    // Определяем объем (цилиндр)
    std::vector<int> all_surfaces = {bottom_surface_tag, top_surface_tag};
    all_surfaces.insert(all_surfaces.end(), side_surface_tags.begin(), side_surface_tags.end());

    int surface_loop_tag = 1;
    gmsh::model::geo::addSurfaceLoop(all_surfaces, surface_loop_tag);
    gmsh::model::geo::addVolume({surface_loop_tag}, 1);



    gmsh::model::geo::synchronize();

    // Добавляем физическую группу (обязательно для сохранения)
    std::vector<int> all_volumes = {1};
    gmsh::model::addPhysicalGroup(3, all_volumes, 1);
    gmsh::model::setPhysicalName(3, 1, "Cylinder");

    gmsh::model::mesh::generate(3);
    gmsh::write("cylinder.msh");

    std::set<std::string> args(argv, argv + argc);
    if (!args.count("-nopopup")) gmsh::fltk::run();

    gmsh::finalize();
    return 0;
}
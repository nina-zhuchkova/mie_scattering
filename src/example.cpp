#include <gmsh.h>
#include <vector>
#include <iostream>

int main() {
    gmsh::initialize();
    gmsh::model::add("scattering");

    // Радиусы
    double R_particle = 1.0;  // Радиус частицы
    double R_outer = 3.0;     // Радиус внешней области
    double lc = 0.01;         // Размер элемента сетки

    // Создаём сферическую частицу
    int particle = gmsh::model::occ::addSphere(0, 0, 0, R_particle);

    // Создаём внешнюю сферическую оболочку
    int outer = gmsh::model::occ::addSphere(0, 0, 0, R_outer);

    // Контейнеры для результатов
    std::vector<std::pair<int, int>> out_dimtags;
    std::vector<std::vector<std::pair<int, int>>> tool_map;

    // Вычитаем частицу из внешней оболочки, чтобы создать окружающую область
    gmsh::model::occ::cut({{3, outer}}, {{3, particle}}, out_dimtags, tool_map, 0, true, true);

    // Синхронизация (чтобы изменения появились в модели)
    gmsh::model::occ::synchronize();

    // Задаём размеры сетки
    gmsh::model::mesh::setSize({{3, outer}}, lc);
    
    // Генерируем 3D-сетку
    gmsh::model::mesh::generate(3);

    // Сохраняем в файл
    gmsh::write("scattering.msh");

    // Отображаем Gmsh GUI (можно закомментировать, если не нужно)
    gmsh::fltk::run();

    // Завершаем работу
    gmsh::finalize();
    return 0;
}

#include <gmsh.h>
#include <iostream>

int main() {
    gmsh::initialize();
    gmsh::model::add("example");

    // Создаём точку
    double x = 0, y = 0, z = 0, cl = 0.1;
    gmsh::model::geo::addPoint(x, y, z, cl, 1);

    // Генерируем mesh
    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate(2);

    gmsh::write("mesh.msh");
    gmsh::finalize();

    std::cout << "Mesh saved to mesh.msh\n";
    return 0;
}

#include <gmsh.h>
#include <iostream>

int main() {
    gmsh::initialize();
    std::cout << "Gmsh успешно инициализирован!" << std::endl;
    gmsh::finalize();
    return 0;
}

cmake_minimum_required(VERSION 3.10)
project(MyGmshProject)

set(CMAKE_CXX_STANDARD 17)

# Указываем путь к Gmsh
set(GMSH_DIR "/usr/local/gmsh-4.11.1-Linux64/lib/cmake/gmsh")
find_package(Gmsh REQUIRED)

add_executable(my_gmsh_app main.cpp)

# Линкуем с Gmsh
target_link_libraries(my_gmsh_app PRIVATE Gmsh::gmsh)

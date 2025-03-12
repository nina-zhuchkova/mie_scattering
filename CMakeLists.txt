cmake_minimum_required(VERSION 3.10)
project(gmsh_only)

# Указываем стандарт C++
set(CMAKE_CXX_STANDARD 11)

# Используем переменную окружения для GMSH
if(DEFINED ENV{GMSH_DIR})
    set(GMSH_DIR $ENV{GMSH_DIR})
endif()
# Находим GMSH
find_package(GMSH REQUIRED)

# Добавляем исходные файлы
add_executable(test_gmsh test_gmsh.cpp)
target_link_libraries(test_gmsh GMSH::GMSH)

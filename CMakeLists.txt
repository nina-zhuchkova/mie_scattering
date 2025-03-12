cmake_minimum_required(VERSION 3.10)
project(MyGmshProject)

set(CMAKE_CXX_STANDARD 17)

# Явно указываем пути к заголовкам и библиотеке
include_directories(/usr/include/gmsh)
link_directories(/usr/lib)

# Добавляем исполняемые файлы
add_executable(main src/main.cpp)
add_executable(example src/example.cpp)

# Линкуем с Gmsh
target_link_libraries(main gmsh)
target_link_libraries(example gmsh)

# Добавляем цель для запуска программ
add_custom_target(run ALL
    COMMAND ./main && ./example
    DEPENDS main example
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
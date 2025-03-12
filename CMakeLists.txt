cmake_minimum_required(VERSION 3.10)
project(MyGmshProject)

set(CMAKE_CXX_STANDARD 17)

# Передаем переменную окружения DISPLAY для всех команд в CMake
set(ENV{DISPLAY} ":0")

# Указываем пути к заголовкам и библиотеке
include_directories(/usr/include/gmsh)
link_directories(/usr/lib)

# Добавляем исполняемые файлы
add_executable(main src/main.cpp)
add_executable(example src/example.cpp)

# Линкуем с GMSH
target_link_libraries(main gmsh)
target_link_libraries(example gmsh)

# Добавляем цель для компиляции и запуска, чтобы сначала запускался main, а затем example
# Использование bash для продолжения работы контейнера
add_custom_target(run
    COMMAND ./build/main && ./build/example
    DEPENDS main example
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)


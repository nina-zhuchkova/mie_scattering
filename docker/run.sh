#!/bin/bash

# Компиляция кода (если нужно)
g++ -o main main.cpp -I/usr/include/gmsh -lgmsh

# Запуск GMSH с файлом model.geo
gmsh -2 model.geo

# Если нужно выполнить что-то еще, добавь сюда
# Например, запуск программы или другие команды

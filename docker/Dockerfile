# Используем официальный образ Ubuntu
FROM ubuntu:22.04

# Устанавливаем зависимости
RUN apt update && apt install -y \
    libglu1-mesa-dev \
    freeglut3-dev \
    mesa-common-dev \
    && rm -rf /var/lib/apt/lists/*

# Устанавливаем Gmsh для C++
RUN apt update && apt install -y gmsh

# Устанавливаем Gmsh для Python
RUN pip3 install gmsh

# Копируем файлы из папки src в контейнер
COPY ./src /app/src

# Устанавливаем точку входа в контейнер
CMD ["/bin/bash"]

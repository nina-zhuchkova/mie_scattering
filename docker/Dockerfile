# Используем официальный образ Ubuntu
FROM ubuntu:22.04

# Устанавливаем зависимости
RUN apt update && apt install -y \
    python3 \
    python3-pip \
    cmake \
    g++ \
    gmsh \
    && rm -rf /var/lib/apt/lists/*

RUN ln -s /usr/bin/python3 /usr/bin/python
# Устанавливаем GMSH для Python
RUN pip3 install gmsh

# Устанавливаем переменную окружения для GMSH
ENV GMSH_DIR=/usr

# Проверка версий CMake и GMSH
RUN cmake --version
RUN gmsh --version

# Копируем файлы из папки src в контейнер
COPY ./src /app/src

# Копируем CMakeLists.txt
COPY ./CMakeLists.txt /app/

# Переходим в директорию с проектом и выполняем сборку
WORKDIR /app

# Создаем директорию для сборки и выполняем сборку в одном RUN шаге
RUN mkdir -p build && cd build && cmake .. && make

# Устанавливаем точку входа в контейнер
CMD ["/bin/bash"]
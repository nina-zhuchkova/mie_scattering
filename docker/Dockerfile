FROM ubuntu:latest

# Устанавливаем необходимые пакеты для компиляции, gmsh и другие зависимости
RUN apt-get update && apt-get install -y \
    x11-apps \
    g++ \
    make \
    cmake \
    gmsh \
    libgmsh-dev \     
    && rm -rf /var/lib/apt/lists/*

# Создаём рабочую директорию
WORKDIR /app

# Копируем CMakeLists.txt из корня проекта
COPY CMakeLists.txt /app/

# Копируем исходники из src/
COPY ./src /app/src

# Создаём папку для сборки
RUN mkdir build
WORKDIR /app/build

# Запускаем CMake и сборку
RUN cmake .. && make

ENV DISPLAY=:0

RUN apt-get install -y x11-xserver-utils

# Запускаем программу
CMD ["./main"]

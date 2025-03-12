FROM ubuntu:latest

# Устанавливаем необходимые пакеты для компиляции, gmsh и другие зависимости
RUN apt-get update && apt-get install -y \
    x11-apps \
    libx11-dev \
    libfltk1.3 \
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
COPY ./docker/run.sh /run.sh

# Создаём папку для сборки
RUN mkdir build
WORKDIR /app/build

ENV DISPLAY=:0

# Запускаем CMake и сборку
RUN cmake .. && make


RUN chmod +x /run.sh
CMD ["/run.sh"]


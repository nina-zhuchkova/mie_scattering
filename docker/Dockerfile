# Используем официальный образ Ubuntu 22.04
FROM ubuntu:22.04  

# Установка зависимостей
RUN apt update && apt install -y \
    build-essential cmake gdb \
    mesa-utils x11-apps libx11-6 \
    git wget vim  

# Устанавливаем Gmsh вручную
RUN wget -q https://gmsh.info/bin/Linux/gmsh-4.11.1-Linux64.tgz && \
    tar -xzf gmsh-4.11.1-Linux64.tgz -C /usr/local && \
    rm gmsh-4.11.1-Linux64.tgz

# Устанавливаем переменные среды для Gmsh
ENV PATH="/usr/local/gmsh-4.11.1-Linux64/bin:$PATH"
ENV GMSH_DIR="/usr/local/gmsh-4.11.1-Linux64/lib/cmake/gmsh"
ENV CMAKE_PREFIX_PATH="/usr/local/gmsh-4.11.1-Linux64/lib/cmake:$CMAKE_PREFIX_PATH"

# Устанавливаем рабочую директорию
WORKDIR /app  

# Копируем исходный код проекта внутрь контейнера
COPY src /app  

# Создаём папку для сборки
RUN mkdir -p build  

# Компиляция проекта
WORKDIR /app/build
RUN cmake .. && make  

# Команда по умолчанию (интерактивная оболочка)
CMD ["bash"]

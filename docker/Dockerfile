FROM ubuntu:latest

# Устанавливаем необходимые пакеты для компиляции, gmsh и другие зависимости
RUN apt-get update && apt-get install -y \
    x11-apps \
    g++ \
    make \
    gmsh \
    libgmsh-dev \     
    && rm -rf /var/lib/apt/lists/*

# Копируем папку src из корневой директории в контейнер
COPY ./src /src
COPY ./docker/run.sh /run.sh

ENV DISPLAY=:0

RUN chmod +x /run.sh

# Переходим в директорию с исходным кодом и компилируем его

CMD ["/run.sh"]

# Запуск gmsh и выполнения программы

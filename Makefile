# Имя образа
IMAGE_NAME = gmsh-env
CONTAINER_NAME = gmsh-container

# Собрать образ
build:
	docker build --no-cache -t $(IMAGE_NAME) -f docker/Dockerfile .

# Запустить контейнер
run:
	docker run --rm -it \
		-e DISPLAY

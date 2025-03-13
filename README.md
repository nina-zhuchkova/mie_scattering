before run

xhost local:docker

to run

sudo docker compose -f ./docker/docker-compose.yml up --build

to interfire at container work

docker exec -ti docker-cpp-gmsh-1 bash

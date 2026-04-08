# mie_scattering

Численный микропроект по `2D` рассеянию на круговой частице на стеке
`gmsh + DOLFINx + UFL + ffcx`.

Репозиторий сейчас ориентирован на современный `DOLFINx`, который удобнее
всего запускать через `docker`-образ `dolfinx/dolfinx:stable`.

На `Ubuntu 22.04` системный `apt`-пакет `dolfinx` оказывается слишком старым
для этого репозитория, поэтому основной workflow проекта считается docker-based.

Сейчас репозиторий сфокусирован на реалистичном минимуме:

- `1.1` плоская волна в `2D`;
- `1.3` гауссов пучок в `2D`.

## Что уже реализовано

- круговая геометрия частицы и внешней области в `gmsh`;
- автоматическая генерация сетки при сборке;
- кусочно-постоянный показатель преломления `n(x)`;
- `2D` скалярная постановка Гельмгольца;
- отдельные расчеты для плоской волны и гауссова пучка;
- сохранение падающего, рассеянного и полного полей.

## Текущая постановка

Сейчас solver работает не через жесткое условие `u = u_inc` на внешней
границе, а через более естественную для scattering-задачи схему:

- неизвестным является рассеянное поле `u_sc`;
- полное поле восстанавливается как `u_tot = u_inc + u_sc`;
- объемный источник получается из падающего поля и материального контраста;
- на внешней границе стоит однородное Robin-условие
  `∂n u_sc + alpha u_sc = 0`.

Это все еще упрощенная модель:

- задача скалярная, а не полная векторная Максвелловская;
- `PETSc` в доступном образе real-valued, поэтому расчет пока тоже
  real-valued;
- Robin-условие здесь играет роль мягкого поглощения на внешней границе,
  а не точного complex-valued radiation condition.

Подробности вынесены в:

- [current_formulation.md](docs/current_formulation.md)
- [mie_implementation_plan.md](docs/mie_implementation_plan.md)

## Структура

- [CMakeLists.txt](CMakeLists.txt):
  сборка проекта, вызов `ffcx`, генерация сетки
- [main.cpp](src/main.cpp):
  чтение сетки, задание коэффициентов, сборка и решение
- [scattering_2d.py](ufl/scattering_2d.py):
  слабая форма `2D` задачи
- [circle_2d.geo](meshes/circle_2d.geo):
  геометрия частицы и внешней области
- [generate_circle_mesh.py](scripts/generate_circle_mesh.py):
  конвертация `gmsh -> DOLFINx/XDMF`

## Сборка

Локально:

```bash
cmake -S . -B build
cmake --build build -j
./build/mie_scattering_2d
```

Если локально `DOLFINx` не установлен, самый надежный путь сейчас такой:

```bash
docker run --rm \
  -v "$PWD":/work \
  -w /work \
  dolfinx/dolfinx:stable \
  bash -lc 'cmake -S . -B build && cmake --build build -j && ./build/mie_scattering_2d'
```

Или короче через готовый скрипт:

```bash
bash scripts/build_docker.sh
```

Скрипт использует отдельный каталог `build-docker`, чтобы не конфликтовать с
локальной сборкой, и запускает контейнер от твоего пользователя.

## Что появляется в results

Сетка:

- `results/circle_2d.xdmf`
- `results/circle_2d.h5`

Материал:

- `results/material_n2_2d.pvd`

Плоская волна:

- `results/plane_wave_2d_incident.pvd`
- `results/plane_wave_2d_scattered.pvd`
- `results/plane_wave_2d_total.pvd`
- `results/plane_wave_2d_background_source.pvd`

Гауссов пучок:

- `results/gaussian_beam_2d_incident.pvd`
- `results/gaussian_beam_2d_scattered.pvd`
- `results/gaussian_beam_2d_total.pvd`
- `results/gaussian_beam_2d_background_source.pvd`

## Что логично делать дальше

- добавить модуль и параметры запуска в командную строку;
- заменить real-valued Robin на более физичную complex-valued постановку;
- сделать сравнение с аналитикой для `2D` цилиндра;
- только потом переносить задачу в `3D`.

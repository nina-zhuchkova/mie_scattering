cmake_minimum_required(VERSION 3.21)

set(PROJECT_NAME mie_scattering)
project(${PROJECT_NAME} LANGUAGES C CXX)

if(NOT TARGET dolfinx)
  find_package(DOLFINX REQUIRED)
endif()

include(CheckSymbolExists)
set(CMAKE_REQUIRED_INCLUDES ${PETSC_INCLUDE_DIRS})
check_symbol_exists(PETSC_USE_COMPLEX petscsystypes.h PETSC_SCALAR_COMPLEX)
check_symbol_exists(PETSC_USE_REAL_DOUBLE petscsystypes.h PETSC_REAL_DOUBLE)

if(PETSC_SCALAR_COMPLEX EQUAL 1)
  if(PETSC_REAL_DOUBLE EQUAL 1)
    set(SCALAR_TYPE "--scalar_type=complex128")
  else()
    set(SCALAR_TYPE "--scalar_type=complex64")
  endif()
else()
  if(PETSC_REAL_DOUBLE EQUAL 1)
    set(SCALAR_TYPE "--scalar_type=float64")
  else()
    set(SCALAR_TYPE "--scalar_type=float32")
  endif()
endif()

add_custom_command(
  OUTPUT scattering_2d.c
  COMMAND ffcx ${CMAKE_CURRENT_SOURCE_DIR}/ufl/scattering_2d.py ${SCALAR_TYPE}
  VERBATIM
  DEPENDS ufl/scattering_2d.py
  COMMENT "Compile scattering_2d.py using FFCx"
)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(MIE_MESH_XDMF ${CMAKE_CURRENT_SOURCE_DIR}/results/circle_2d.xdmf)
set(MIE_MESH_H5 ${CMAKE_CURRENT_SOURCE_DIR}/results/circle_2d.h5)
set(MIE_RESULTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/results)

add_custom_command(
  OUTPUT ${MIE_MESH_XDMF} ${MIE_MESH_H5}
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/results
  COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/scripts/generate_circle_mesh.py
          --geo ${CMAKE_CURRENT_SOURCE_DIR}/meshes/circle_2d.geo
          --xdmf ${MIE_MESH_XDMF}
          --name circle_2d
  DEPENDS meshes/circle_2d.geo scripts/generate_circle_mesh.py
  COMMENT "Generate 2D circle mesh from gmsh geometry"
  VERBATIM
)

add_custom_target(generate_circle_mesh DEPENDS ${MIE_MESH_XDMF} ${MIE_MESH_H5})

add_executable(${PROJECT_NAME}_2d
  src/main.cpp
  ${CMAKE_CURRENT_BINARY_DIR}/scattering_2d.c
)

target_link_libraries(${PROJECT_NAME}_2d dolfinx)
target_compile_features(${PROJECT_NAME}_2d PUBLIC cxx_std_20)
set(CMAKE_CXX_EXTENSIONS OFF)
add_dependencies(${PROJECT_NAME}_2d generate_circle_mesh)
target_compile_definitions(${PROJECT_NAME}_2d PRIVATE
  MIE_MESH_XDMF="${MIE_MESH_XDMF}"
  MIE_RESULTS_DIR="${MIE_RESULTS_DIR}"
)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-Wno-comment" HAVE_NO_MULTLINE)
set_source_files_properties(
  src/main.cpp
  PROPERTIES
    COMPILE_FLAGS
    "$<$<BOOL:${HAVE_NO_MULTLINE}>:-Wno-comment -Wall -Wextra -pedantic>"
)

#!/usr/bin/env python3
"""Generate a 2D circle mesh from a gmsh .geo file and store it as XDMF."""

from __future__ import annotations

import argparse
from pathlib import Path

import gmsh  # type: ignore
from mpi4py import MPI

from dolfinx import io
from dolfinx.io import XDMFFile


def write_mesh(name: str, xdmf_path: Path, mesh_data) -> None:
    mesh_data.mesh.name = name
    if mesh_data.cell_tags is not None:
        mesh_data.cell_tags.name = f"{name}_cells"
    if mesh_data.facet_tags is not None:
        mesh_data.facet_tags.name = f"{name}_facets"

    tdim = mesh_data.mesh.topology.dim
    for dim in range(tdim):
        mesh_data.mesh.topology.create_connectivity(dim, tdim)

    geometry_xpath = f"/Xdmf/Domain/Grid[@Name='{name}']/Geometry"
    with XDMFFile(mesh_data.mesh.comm, str(xdmf_path), "w") as file:
        file.write_mesh(mesh_data.mesh)
        if mesh_data.cell_tags is not None:
            file.write_meshtags(
                mesh_data.cell_tags,
                mesh_data.mesh.geometry,
                geometry_xpath=geometry_xpath,
            )
        if mesh_data.facet_tags is not None:
            file.write_meshtags(
                mesh_data.facet_tags,
                mesh_data.mesh.geometry,
                geometry_xpath=geometry_xpath,
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--geo", required=True, help="Path to gmsh .geo file")
    parser.add_argument("--xdmf", required=True, help="Output XDMF path")
    parser.add_argument("--name", default="circle_2d", help="Mesh name in XDMF")
    parser.add_argument("--dim", type=int, default=2, help="Topological mesh dimension")
    parser.add_argument("--gdim", type=int, default=2, help="Geometric mesh dimension")
    args = parser.parse_args()

    geo_path = Path(args.geo).resolve()
    xdmf_path = Path(args.xdmf).resolve()
    xdmf_path.parent.mkdir(parents=True, exist_ok=True)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.open(str(geo_path))
    gmsh.model.mesh.generate(args.dim)
    gmsh.model.mesh.setOrder(1)

    mesh_data = io.gmsh.model_to_mesh(gmsh.model, MPI.COMM_WORLD, 0, gdim=args.gdim)
    write_mesh(args.name, xdmf_path, mesh_data)
    gmsh.finalize()


if __name__ == "__main__":
    main()

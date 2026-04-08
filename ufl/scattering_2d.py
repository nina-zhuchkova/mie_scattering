"""Scalar 2D scattered-field form for the mie_scattering project."""

from basix.ufl import element
from ufl import (
    Coefficient,
    Constant,
    FunctionSpace,
    Mesh,
    TestFunction,
    TrialFunction,
    ds,
    dx,
    grad,
    inner,
)

element_u = element("Lagrange", "triangle", 1)
element_q = element("DG", "triangle", 0)
coord_element = element("Lagrange", "triangle", 1, shape=(2,))

mesh = Mesh(coord_element)
V = FunctionSpace(mesh, element_u)
Q = FunctionSpace(mesh, element_q)

u = TrialFunction(V)
v = TestFunction(V)

n2 = Coefficient(Q)
u_inc = Coefficient(V)
source_bg = Coefficient(V)

k0_sq = Constant(mesh)
n_background_sq = Constant(mesh)
alpha = Constant(mesh)

a = inner(grad(u), grad(v)) * dx - k0_sq * n2 * inner(u, v) * dx + alpha * inner(u, v) * ds
L = inner(source_bg, v) * dx + k0_sq * (n2 - n_background_sq) * inner(u_inc, v) * dx

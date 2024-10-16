import dolfin
import numpy as np
from fenics import Mesh, TestFunctions, Constant, assemble, Expression, interpolate, solve, DirichletBC, plot, errornorm, set_log_active, derivative, parameters, split, dot, div, CompiledSubDomain, MeshFunction, sqrt, Measure, FacetNormal, Identity
from solver import *
from laser.laser_source_mesher import save_laser_mesh

vfile = File("laser/data/laser_v.pvd", "compressed")
efile = File("laser/data/laser_e.pvd", "compressed")

parameters = {
    "s_v": 2,
    "s_p": 1,
    "s_e": 1,
    "r": 1,
    "start_time": 0.,
    "end_time":  0.4,
    "slab_size": 0.001,
    "n_x": 1,
    "nu0": 1.0,
    "alpha": 6.88e-5,
    "g": 9.81,
    "R": 1.0,
    "T0": 200,
    "E_A": 1.0,
    "rho": 1.0,
    "k": 1.0,
}

# get spatial function space
space_mesh = Mesh("laser/fine_laser.xml")

# boundaries
walls = CompiledSubDomain("(near(x[0], 0) || near(x[1], 0) || near(x[0], 1) || near(x[1], 1)) && on_boundary")

facet_marker = MeshFunction("size_t", space_mesh, 1)
facet_marker.set_all(0)
walls.mark(facet_marker, 1)

# define time dependent Dirichlet boundary conditions
def get_bcs(V, Time):
    bcs = []
    offset = 2*len(Time.dof_locations)
    for i, t_q in enumerate(Time.dof_locations):
        bcs.append(DirichletBC(V.sub(i), Constant((0, 0)), walls))
        bcs.append(DirichletBC(V.sub(i + offset), Constant(293.15), walls))
    return bcs

def avg_temp(Time, v0, p0, e0):
    dx = Measure('dx', domain=space_mesh)
    return assemble(e0*dx)

x0 = (0.75, 0.75)
laser = Expression(f'pow(10,5)*sqrt(2*pi)*exp(-pow(10, 4) * (pow(x[0]-0.75, 2) + pow(x[1] - 0.75, 2)))', degree=2)

sim = Boussinesque_Solver('laser', space_mesh, parameters, inhomogeneity = laser)
sim.set_goal_functional(avg_temp)
sim.solve(get_bcs, vfile, efile, initial_condition=Constant((0.,0.,0.,0.)), QoI = avg_temp)
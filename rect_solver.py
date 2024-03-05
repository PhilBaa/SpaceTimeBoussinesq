import dolfin
import numpy as np
from fenics import Mesh, TestFunctions, Constant, assemble, Expression, interpolate, solve, DirichletBC, plot, errornorm, set_log_active, derivative, parameters, split, dot, div, CompiledSubDomain, MeshFunction, sqrt, Measure, FacetNormal, Identity
from solver import *
from rect.rect_mesher import save_rect_mesh

vfile = File("rect/data/rect_v.pvd", "compressed")
efile = File("rect/data/rect_e.pvd", "compressed")
start_time = 0.
end_time = 0.4

s_v = 2
s_p = 1
s_e = 2
r = 1
slab_size = 0.002
n_x = 1
nu = 0.001
alpha = 6.88e-5
g = 9.81

parameters = {
    "s_v": s_v,
    "s_p": s_p,
    "s_e": s_e,
    "r": r,
    "start_time": start_time,
    "end_time": end_time,
    "slab_size": slab_size,
    "n_x": n_x,
    "nu": nu,
    "alpha": alpha,
    "g": g
}

space_mesh = Mesh("rect/fine_rect.xml") 

# boundaries
tol = 1e-6

walls = CompiledSubDomain(f"((near(x[0], -2.0) && x[1]>0.5 + {tol}) || near(x[1], 0.0) || near(x[1], 3.0) || (near(x[0], 8.0 && x[1] < 2.5 - {tol}))) && on_boundary")

inflow = CompiledSubDomain(f"near(x[0], -2.0) && x[1]<0.5 + {tol} && on_boundary")
outflow = CompiledSubDomain(f"near(x[0], 8.0) && x[1]>2.5 - {tol} && on_boundary")

rects = CompiledSubDomain(f"((near(x[0], 1.0) && x[1] < 1.0) || \
                            (near(x[0], 2.0) && x[1] < 1.0) || \
                            (near(x[0], 3.0) && x[1] < 2.0) || \
                            (near(x[0], 4.0) && x[1] < 2.0) || \
                            (near(x[1], 1.0) && x[0] > 1.0-{tol} && x[0] < 2.0 + {tol}) || \
                            (near(x[1], 2.0) && x[0] > 3.0-{tol} && x[0] < 4.0 + {tol}) ) \
                            && on_boundary")

def get_bcs(V, Time):
    bcs = []
    offset = 2*len(Time.dof_locations)
    for i, t_q in enumerate(Time.dof_locations):
        bcs.append(DirichletBC(V.sub(i), Constant((0, 0)), walls))
        bcs.append(DirichletBC(V.sub(i), Constant((0, 0)), rects))
        bcs.append(DirichletBC(V.sub(i), Constant((1, 0)), inflow))
        bcs.append(DirichletBC(V.sub(i), Constant((1, 0)), outflow))

        bcs.append(DirichletBC(V.sub(i + offset), Constant(1), walls))
        bcs.append(DirichletBC(V.sub(i + offset), Constant(1), inflow))
        bcs.append(DirichletBC(V.sub(i + offset), Constant(200), rects))
    return bcs

facet_marker = MeshFunction("size_t", space_mesh, 1)
facet_marker.set_all(0)
walls.mark(facet_marker, 1)

sim = Boussinesque_Solver('rect', space_mesh, parameters)
sim.solve(get_bcs, vfile, efile)

vfile = File("rect/data/forget.pvd", "compressed")
efile = File("rect/data/forget.pvd", "compressed")
for i in np.linspace(0, 1, 5):
    save_rect_mesh('rect/meshes/rect' + str(i), i)
    sim = Boussinesque_Solver('rect', Mesh('rect/meshes/rect' + str(i) + '.xml'), parameters)
    sim.solve(get_bcs, vfile, efile)
    print('Done with ' + str(i))


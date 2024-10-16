import dolfin
import numpy as np
from fenics import Mesh, TestFunctions, Constant, assemble, Expression, interpolate, solve, DirichletBC, plot, errornorm, set_log_active, derivative, parameters, split, dot, div, CompiledSubDomain, MeshFunction, sqrt, Measure, FacetNormal, Identity
from solver import *
from rect.rect_mesher import save_rect_mesh
import pandas as pd

vfile = File("rect/data/rect_v2.pvd", "compressed")
efile = File("rect/data/rect_e2.pvd", "compressed")
start_time = 0.
end_time = 40.0

s_v = 2
s_p = 1
s_e = 1
r = 1
slab_size = 0.3
n_x = 1

Re = 50
Pe = 36
Ri = 0

parameters = {
    "s_v": s_v,
    "s_p": s_p,
    "s_e": s_e,
    "r": r,
    "start_time": start_time,
    "end_time": end_time,
    "slab_size": slab_size,
    "n_x": n_x,
    "nu0": Re,
    "alpha": 1,
    "g": Ri,
    "R": 1,
    "T0":1000,
    "E_A": 0,
    "rho": 1/Re,
    "k": 1/Pe,
}

space_mesh = Mesh("rect/fine_rect.xml") 

# boundaries
tol = 1e-6

walls = CompiledSubDomain(f"((near(x[0], -2.0) && x[1]>0.5 + {tol}) || near(x[1], 0.0) || near(x[1], 3.0) || (near(x[0], 8.0) && x[1] < 2.5 - {tol})) && on_boundary")

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

        bcs.append(DirichletBC(V.sub(i + offset), Constant(0), walls))
        bcs.append(DirichletBC(V.sub(i + offset), Constant(0), inflow))
        bcs.append(DirichletBC(V.sub(i + offset), Constant(1), rects))
    return bcs

def Temperature_at_outflow(Time, solutions_v, solutions_p, e0):
    facet_marker = MeshFunction("size_t", space_mesh, 1)
    facet_marker.set_all(0)
    outflow.mark(facet_marker, 2)
    ds = Measure("ds", domain=space_mesh, subdomain_data=facet_marker)
    return assemble(e0*ds(2))

facet_marker = MeshFunction("size_t", space_mesh, 1)
facet_marker.set_all(0)
walls.mark(facet_marker, 1)

sim = Boussinesque_Solver('rect', space_mesh, parameters)
#sim.set_lagrange_multiplier(1000)
sim.set_goal_functional(Temperature_at_outflow)

sim.solve(get_bcs, vfile, efile, initial_condition=Constant((0.0,0.0,0.0,0.0)))

t = np.arange(start_time, end_time, 0.3)

#np.savetxt("rect/data/outflow_temp_1.csv", sim.func_vals)

vfile = File("rect/data/forget.pvd", "compressed")
efile = File("rect/data/forget.pvd", "compressed")

res = np.linspace(0.1, 0.9, 10)

data = np.zeros((len(res), len(t)))
dofs = []

for i, r in enumerate(res):
    save_rect_mesh('rect/meshes/rect' + str(i), r)
    sim = Boussinesque_Solver('rect', Mesh('rect/meshes/rect' + str(i) + '.xml'), parameters)
    sim.set_lagrange_multiplier(1000)
    sim.set_goal_functional(Temperature_at_outflow)
    try:
        sim.solve(get_bcs, vfile, efile, initial_condition=Constant((0.0,0.0,0.0,0.0)))
        data[i] = sim.func_vals
    except RuntimeError:
        pass
    finally:
        dofs.append(sim.Vh.dim())

pd.DataFrame(data, index=dofs, columns = t).to_csv("rect/data/outflow_temp_dofs.csv")



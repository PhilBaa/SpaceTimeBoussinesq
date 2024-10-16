import dolfin
import numpy as np
from fenics import Mesh, TestFunctions, Constant, assemble, Expression, interpolate, solve, DirichletBC, plot, errornorm, set_log_active, derivative, parameters, split, dot, div, CompiledSubDomain, MeshFunction, sqrt, Measure, FacetNormal, Identity
from solver import *

vfile = File("schaefer_turek/data/schaefer_turek_v.pvd", "compressed")
pfile = File("schaefer_turek/data/schaefer_turek_p.pvd", "compressed")
efile = File("schaefer_turek/data/schaefer_turek_e.pvd", "compressed")

parameters = {
    "s_v": 2,
    "s_p": 1,
    "s_e": 1,
    "r": 1,
    "start_time": 0.,
    "end_time":  8.0,
    "slab_size": 0.01,
    "n_x": 1,
    "nu0": 0.001,
    "alpha": 5.0e-3,
    "g": 9.81,
    "R": 1.0,
    "T0": 200,
    "E_A": 1.0 ,
    "rho": 1.0,
    "k": 0.005,
}

# get spatial function space
space_mesh = Mesh("schaefer_turek/schaefer_turek.xml")

# boundaries
inflow = CompiledSubDomain("near(x[0], 0) && on_boundary")
outflow = CompiledSubDomain("near(x[0], 2.2) && on_boundary")
walls = CompiledSubDomain("near(x[1], 0) || near(x[1], 0.41) && on_boundary")
cylinder = CompiledSubDomain("x[0]>0.1 && x[0]<0.3 && x[1]>0.1 && x[1]<0.3 && on_boundary")

facet_marker = MeshFunction("size_t", space_mesh, 1)
facet_marker.set_all(0)
inflow.mark(facet_marker, 1)
outflow.mark(facet_marker, 2)
walls.mark(facet_marker, 3)
cylinder.mark(facet_marker, 4)

ds_cylinder = Measure("ds", subdomain_data=facet_marker, subdomain_id=4)

# define time dependent Dirichlet boundary conditions
def get_bcs(V, Time):
    bcs = []
    offset = 2*len(Time.dof_locations)
    for i, t_q in enumerate(Time.dof_locations):
        inflow_parabola = ("5.0*sin(0.125*pi*t)*x[1]*(0.41 - x[1]) / pow(0.41, 2)", "0")
        bcs.append(DirichletBC(V.sub(i), Expression(inflow_parabola, degree=2, pi=np.pi, t=t_q), inflow))
        bcs.append(DirichletBC(V.sub(i), Constant((0, 0)), walls))
        bcs.append(DirichletBC(V.sub(i), Constant((0, 0)), cylinder))
        bcs.append(DirichletBC(V.sub(i + offset), Constant(600), cylinder))
        bcs.append(DirichletBC(V.sub(i + offset), Constant(0), inflow))
        bcs.append(DirichletBC(V.sub(i + offset), Constant(0), walls))
    return bcs

sim = Boussinesque_Solver('schaefer_turek', space_mesh, parameters)





def drag_and_lift(self, Time, v0, p0, e0):
    dU = TrialFunction(self.Vh) 
    dv, dp, de = split(dU)
    nu = self.nu0 * exp(self.E_A/(self.R*(e0 + self.T0))) 
    D = 0.1
    v_bar = 2/3 * 4.0*1.5*0.205*(0.41 - 0.205) / pow(0.41, 2)

    n = FacetNormal(space_mesh)

    drag = assemble(
        2/(v_bar**2*D)*
        (
        - dot(p0 * Identity(len(v0)) , -n)[0]
        + nu * dot(grad(v0), -n)[0]
        ) * ds_cylinder
    )

    lift = assemble(
        2/(v_bar**2*D)*
        (
        - dot(p0 * Identity(len(v0)) , -n)[1]
        + nu * dot(grad(v0), -n)[1]
        ) * ds_cylinder
    )
    return drag, lift

sim.set_QoI(drag_and_lift)
sim.solve(get_bcs, vfile, pfile, efile, initial_condition=Constant((0.,0.,0.,0.)))

drag_lift_values = np.asarray(sim.get_QoI_values())

drag = drag_lift_values[:,0]
lift = drag_lift_values[:,1]

np.savetxt('schaefer_turek/data/drag_lift_values.txt', np.column_stack((drag, lift)))


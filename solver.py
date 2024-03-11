import dolfin
import numpy as np
#import matplotlib.pyplot as plt
from sympy import Mul, Id, symbols, init_printing, expand, compose, diff, lambdify, Piecewise
from IPython.display import Math, display
from fenics import Mesh, VectorElement, Function, TrialFunction, TestFunction, TestFunctions, FunctionSpace, dx, inner, grad, FiniteElement, MixedElement, Constant, assemble, Expression, interpolate, solve, DirichletBC, plot, errornorm, set_log_active, derivative, parameters, split, dot, div, CompiledSubDomain, MeshFunction, sqrt, Measure, FacetNormal, Identity
from ufl import replace
from ufl.operators import exp
from fenics import exp
import time
import pandas as pd

from dolfin.cpp.io import File
set_log_active(False) # turn off FEniCS logging
parameters["reorder_dofs_serial"] = False # to avoid reordering of dofs
t = symbols("t")
init_printing()

TIME_DG_QUADRATURE = "Gauss-Lobatto"

# HELPER FUNCTIONS
# compute temporal basis functions from the roots
def compute_basis_functions(roots):
    basis = []
    for i in range(len(roots)):
        f = 1 + 0*Id(t)
        for j in range(len(roots)):
            if i != j:
                f = Mul(f,(t - roots[j]) / (roots[i] - roots[j]))
        #basis.append(f)
        basis.append(expand(f))
    return basis

# transform roots from [-1, 1] to [0, 1]
def transform_roots(roots):
    new_roots = []
    for root in roots:
        new_roots.append((root + 1) / 2)
    return new_roots

# dictionary of temporal bases depending on the polynomial degree
# FE basis on reference element (0,1)
φ = {}
roots = {}

roots[0] = [1.]
φ[0] = [1 + 0*Id(t)]
roots[1] = [0., 1.]
φ[1] = [1-t, t]
roots[2] = transform_roots([-1, 0, 1])
φ[2] = compute_basis_functions(roots[2])

class TimeFE:
    def __init__(self, r=1, a=0., b=1., n_time=1, n_q_points=3):
        self.r = r # polynomial degree in time
        self.n_time = n_time # number of temporal elements
        self.n_dofs = (self.r+1) * self.n_time # total number of temporal dofs on slab
        self.dof_locations = []
        self.a = a # start time of slab
        self.b = b # end time of slab
        self.k = (self.b - self.a) / self.n_time
        self.epsilon = self.k * 1e-10

        self.generate_mesh()
        self.get_full_basis()
        self.lambdify_basis()
        self.compute_quadrature(n_q_points)

    def generate_mesh(self):
        # create an uniform temporal mesh with mesh size self.k
        self.mesh = [(self.a, self.a+self.k)]
        while len(self.mesh) < self.n_time:
            self.mesh.append((self.mesh[-1][1], self.mesh[-1][1]+self.k))

    # transform basis functions from [0,1] to [a,b]
    def transform_function(self, f, a, b):
        return compose(f, (t-a)/(b-a)) # = f((t-a)/(b-a))
    
    def transform_derivative(self, a, b):
        return 1 / (b-a)

    # get full FE basis and its derivative on temporal mesh
    def get_full_basis(self):
        self._basis = []
        self._basis_derivative = []
        self.local_dofs = {}
        i = 0
        for (a,b) in self.mesh:
            self.local_dofs[(a,b)] = []
            for f, t_q in zip(φ[self.r], roots[self.r]):
                self._basis.append(self.transform_function(f, a, b))
                self._basis_derivative.append(diff(self._basis[-1],t))
                #print(diff(self._basis[-1],t))
                #print(self.transform_function(diff(f,t), a, b) * self.transform_derivative(a, b))
                self.local_dofs[(a,b)].append(i)
                self.dof_locations.append(t_q*(b-a)+a)
                i += 1

    # convert SymPy functions to Python functions and ensure that they are 0 outside the element that they are defined on
    def lambdify_basis(self):
        self.phi = []
        self.dt_phi = []

        for (a,b) in self.mesh:
            for i in self.local_dofs[(a,b)]:
                self.phi.append(
                    lambdify(
                        t,
                        Piecewise(
                            (0, t < a),
                            (0, t > b),
                            (self._basis[i], True)
                        )
                    )
                )

                self.dt_phi.append(
                    lambdify(
                        t,
                        Piecewise(
                            (0, t < a),
                            (0, t > b),
                            (self._basis_derivative[i], True)
                        )
                    )
                )

    def compute_quadrature(self, n_q_points):
        # Gauss-Legendre quadrature points and weights on [-1,1]
        quad_points, quad_weights = np.polynomial.legendre.leggauss(n_q_points)

        # transform quadrature points and weights from [-1,1] to [a,b] for each temporal element [a,b]
        self.quadrature = {}
        for (a, b) in self.mesh:
            t_q = 0.5 * (b-a) * quad_points + 0.5 * (a+b)
            w_q = 0.5 * (b-a) * quad_weights
            self.quadrature[(a,b)] = [(t_q[i], w_q[i]) for i in range(t_q.shape[0])]
            
        # Gauss-Legendre quadrature points and weights on [-1,1]
        quad_points, quad_weights = np.polynomial.legendre.leggauss(n_q_points+2)

        # transform quadrature points and weights from [-1,1] to [a,b] for each temporal element [a,b]
        self.quadrature_fine = {}
        for (a, b) in self.mesh:
            t_q = 0.5 * (b-a) * quad_points + 0.5 * (a+b)
            w_q = 0.5 * (b-a) * quad_weights
            self.quadrature_fine[(a,b)] = [(t_q[i], w_q[i]) for i in range(t_q.shape[0])]

    def plot_basis(self, basis_type="function", title=None):
        assert basis_type in ["function", "derivative"], f"basis_type='{basis_type}' has not been implemented."

        _t = np.linspace(self.a, self.b, 100)
        for i in range(len(self.phi)):
            if basis_type == "function":
                plt.plot(_t, [self.phi[i](time_point) for time_point in _t], label=rf"$\varphi_{{{i}}}$")
            elif basis_type == "derivative":
                plt.plot(_t, [self.dt_phi[i](time_point) for time_point in _t], label=rf"$\partial_t \varphi_{{{i}}}$")

        plt.xlabel("t")
        plt.ylabel("f(t)")
        plt.legend()
        if title is not None:
            plt.title(title)
        plt.show()
        
    def get_solution_at_time(self, time_point, U):
        tmp = np.zeros_like(U[0])
        for i in range(self.n_dofs):
            tmp += self.phi[i](time_point) * U[i]
        return tmp


class Boussinesque_Solver:
    def __init__(self, name, space_mesh, parameters, inhomogeneity = None):
        self.name = name
        self.space_mesh = space_mesh
        self.inhomogeneity = inhomogeneity
        self.read_parameters(parameters)

        self.slabs = [(self.start_time, self.start_time+self.slab_size)]
        while self.slabs[-1][1] < self.end_time - 1e-8:
            self.slabs.append((self.slabs[-1][1], self.slabs[-1][1]+self.slab_size))

        self.set_function_spaces()
        self.goal_functional = None
        self.func_vals = []

        self.lme = 1
    

    def read_parameters(self, parameters):
        self.parameters = parameters
        self.nu0 = parameters["nu0"]
        self.alpha = parameters["alpha"]
        self.g = parameters["g"]
        self.R = parameters["R"]
        self.E_A = parameters["E_A"]
        self.T0 = parameters["T0"]
        self.k = parameters["k"]
        self.rho = parameters["rho"]


        self.s_v = parameters["s_v"]
        self.s_p = parameters["s_p"]
        self.s_e = parameters["s_e"]
        self.r = parameters["r"]
        self.start_time = parameters["start_time"]
        self.end_time = parameters["end_time"]
        self.slab_size = parameters["slab_size"]

    
    def set_lagrange_multiplier(self, lme):
        self.lme = lme


    def set_function_spaces(self):
        self.element = {
            "v": VectorElement("Lagrange", self.space_mesh.ufl_cell(), self.s_v),
            "p": FiniteElement("Lagrange", self.space_mesh.ufl_cell(), self.s_p),
            "e": FiniteElement("Lagrange", self.space_mesh.ufl_cell(), self.s_e),
        }
        Vh = FunctionSpace(self.space_mesh, MixedElement(*self.element.values())) # spatial function space for a single time point
        print(f"Number of spatial DoFs: {Vh.dim()} ({Vh.sub(0).dim()} + {Vh.sub(1).dim()} + {Vh.sub(2).dim()})")
        self.Uh = Function(Vh)
        self.Phih = TestFunctions(Vh)
        self.Vh = Vh

    
    def set_goal_functional(self, goal_functional):
        self.goal_functional = goal_functional
    

    def get_goal_functional_values(self):
        return self.func_vals

    
    def weak_formulation(self, U, Phi, Time, v0, p0, e0, k):
        # start with "empty" space-time variational form
        F = Constant(0.)*U["p"][0]*Phi["p"][0]*dx

        # ================= #
        #   (v,v) - Block   #
        # ================= #

        # volume integrals
        for time_element in Time.mesh:
            # assemble linear terms
            for i in Time.local_dofs[time_element]:
                for j in Time.local_dofs[time_element]:
                    for (t_q, w_q) in Time.quadrature[time_element]:
                        # TODO: to reduce the number of terms in the sum, the sum over the temporal quadrature can be evaluated prior to adding to the form F
                        F += Constant(self.rho * w_q * Time.dt_phi[j](t_q) * Time.phi[i](t_q)) \
                            * dot(U["v"][j], Phi["v"][i]) * dx
            # assemble nonlinearity
            for i in Time.local_dofs[time_element]:
                for j in Time.local_dofs[time_element]:
                    for l in Time.local_dofs[time_element]:
                        # NOTE: For nonlinearities make sure that the temporal quadrature is fine enough
                        # E.g. for the nonlinearity u^2, we need to be able to integrate polynomials of degree 3r exactly in time
                        #      for this we need Gauss-Legendre quadrature of degree >= (3r+1)/2
                        for (t_q, w_q) in Time.quadrature_fine[time_element]:
                            F += Constant(w_q * Time.phi[j](t_q) * Time.phi[l](t_q) * Time.phi[i](t_q)) \
                                *Constant(self.rho)* dot(dot(grad(U["v"][j]), U["v"][l]), Phi["v"][i])* dx   

            for i in Time.local_dofs[time_element]:
                for j in Time.local_dofs[time_element]:
                        for (t_q, w_q) in Time.quadrature_fine[time_element]:
                            theta_tq = 0. # Function?
                            for l in Time.local_dofs[time_element]:
                                theta_tq += exp(self.E_A/(self.R*(U["e"][l]*Time.phi[l](t_q) + self.T0)))
                            F += Constant(self.nu0 * w_q * Time.phi[j](t_q)  * Time.phi[i](t_q)) \
                                * inner(grad(U["v"][j]), grad(Phi["v"][i]))*theta_tq * dx 



        # RHS integral
        for n, time_element in enumerate(Time.mesh):
            for i in Time.local_dofs[time_element]:
                # initial condition
                if k == 0:
                    F -=  Constant(self.lme*Time.phi[i](time_element[0]+Time.epsilon)) * inner(v0, Phi["v"][i]) * dx

        # jump terms (NOTE: For Gauss-Lobatto dG basis, we could hard code the values for the jump term)
        for n, time_element in enumerate(Time.mesh):
            # a) v_m^+ * φ_m^{v,+}
            for i in Time.local_dofs[time_element]:
                for j in Time.local_dofs[time_element]:
                    F += Constant(self.lme*Time.phi[j](time_element[0]+Time.epsilon) * Time.phi[i](time_element[0]+Time.epsilon)) * inner(U["v"][j], Phi["v"][i]) * dx

            # b) v_{m-1}^- * φ_m^{v,+}
            if n > 0:
                prev_time_element = Time.mesh[n-1]
                for i in Time.local_dofs[time_element]:
                    for j in Time.local_dofs[prev_time_element]:
                        F += Constant((-self.lme) * Time.phi[j](prev_time_element[1]-Time.epsilon) * Time.phi[i](time_element[0]+Time.epsilon)) * inner(U["v"][j], Phi["v"][i]) * dx       

        # ================= #
        #   (v,p) - Block   #
        # ================= #

        # volume integral
        for time_element in Time.mesh:
            for i in Time.local_dofs[time_element]:
                for j in Time.local_dofs[time_element]:
                    for (t_q, w_q) in Time.quadrature[time_element]:
                        F += Constant(w_q * Time.phi[j](t_q) * Time.phi[i](t_q)) \
                            * Constant(-1.) * U["p"][j] * div(Phi["v"][i]) * dx

        # ================= #
        #   (p,v) - Block   #
        # ================= #

        # volume integral
        for time_element in Time.mesh:
            for i in Time.local_dofs[time_element]:
                for j in Time.local_dofs[time_element]:
                    for (t_q, w_q) in Time.quadrature[time_element]:
                        F += Constant(w_q * Time.phi[j](t_q) * Time.phi[i](t_q)) \
                            * Constant(-1.) * div(U["v"][j]) * Phi["p"][i] * dx

        # ================= #
        #   (e,e) - Block   #
        # ================= #

        # volume integral
        for time_element in Time.mesh:
            for i in Time.local_dofs[time_element]:
                for j in Time.local_dofs[time_element]:
                    for (t_q, w_q) in Time.quadrature[time_element]:
                        F += Constant(self.k * w_q * Time.phi[j](t_q) * Time.phi[i](t_q)) \
                            * inner(grad(U["e"][j]), grad(Phi["e"][i])) * dx
                        F += Constant(w_q * Time.dt_phi[j](t_q) * Time.phi[i](t_q)) \
                            * U["e"][j] * Phi["e"][i] * dx

        # RHS integral
        for n, time_element in enumerate(Time.mesh):
            print(time_element)
            for i in Time.local_dofs[time_element]:
                # initial condition
                if n == 0:
                    F -= Constant(self.lme * Time.phi[i](time_element[0]+Time.epsilon)) * e0 * Phi["e"][i] * dx
                # inhomogeneity
                if self.inhomogeneity is not None:
                    for (t_q, w_q) in Time.quadrature[time_element]:
                        F -= Constant(w_q * Time.phi[i](t_q)) * self.inhomogeneity * Phi["e"][i] * dx #time integral


        #jump terms
        for n, time_element in enumerate(Time.mesh):
            # a) e_m^+ * φ_m^{e,+}
            for i in Time.local_dofs[time_element]:
                for j in Time.local_dofs[time_element]:
                    F += Constant(self.lme * Time.phi[j](time_element[0]+Time.epsilon) * Time.phi[i](time_element[0]+Time.epsilon)) * U["e"][j] * Phi["e"][i] * dx

            # b) v_{m-1}^- * φ_m^{v,+}
            if n > 0:
                prev_time_element = Time.mesh[n-1]
                for i in Time.local_dofs[time_element]:
                    for j in Time.local_dofs[prev_time_element]:
                        F += Constant((-1.0 * self.lme) * Time.phi[j](prev_time_element[1]-Time.epsilon) * Time.phi[i](time_element[0]+Time.epsilon)) * U["e"][j] * Phi["e"][i] * dx       

        # ================= #
        #   (e,v) - Block   #
        # ================= #

        # volume integral
        for time_element in Time.mesh:
            for i in Time.local_dofs[time_element]:
                for j in Time.local_dofs[time_element]:
                    for l in Time.local_dofs[time_element]:
                        for (t_q, w_q) in Time.quadrature[time_element]:
                            F += Constant(w_q * Time.phi[j](t_q) * Time.phi[i](t_q) * Time.phi[l](t_q)) \
                                * inner(U["v"][j], grad(U["e"][i])) * Phi["e"][l] * dx
            #linear term due to thermal expansion
            for i in Time.local_dofs[time_element]:
                for j in Time.local_dofs[time_element]:
                    for (t_q, w_q) in Time.quadrature[time_element]:
                        F += Constant(w_q * Time.phi[j](t_q) * Time.phi[i](t_q)* self.alpha* self.g) \
                            *U["e"][i] * Phi["v"][j][1] * dx
        return F


    def solve(self, get_bcs: callable, vfile, efile, initial_condition = Constant((0.,0.,0.,1.)), goal_functional = None):
        # initial condition on slab
        U0 = Function(self.Vh)
        U1 = Function(self.Vh)
        v0, p0, e0 = split(U0)
        U0 = interpolate(initial_condition, self.Vh)

        ##############################################
        # Start a time marching / time slabbing loop #
        ##############################################

        for k, slab in enumerate(self.slabs):
            print(f"Solving on slab_{k} = Ω x ({round(slab[0],5)}, {round(slab[1],5)}) ...")

            #########################################
            # Create temporal finite element object #
            #########################################
            Time = TimeFE(r=self.r, a=slab[0], b=slab[1], n_time=1, n_q_points=self.r+1)

            V = FunctionSpace(self.space_mesh, MixedElement(*[self.element["v"] for _ in range(Time.n_dofs)], *[self.element["p"] for _ in range(Time.n_dofs)], *[self.element["e"] for _ in range(Time.n_dofs)]))
            n_subspaces = V.num_sub_spaces()

            U_kh = Function(V)
            Phi_kh = TestFunctions(V)

            # split U_kh and Phi_kh into velocity and pressure parts
            U = {"v": split(U_kh)[:Time.n_dofs], "p": split(U_kh)[Time.n_dofs: 2*Time.n_dofs], "e": split(U_kh)[2*Time.n_dofs:]}
            Phi = {"v": Phi_kh[:Time.n_dofs], "p": Phi_kh[Time.n_dofs: 2*Time.n_dofs], "e": Phi_kh[2*Time.n_dofs:]}

            F = self.weak_formulation(U, Phi, Time, v0, p0, e0, k)
            bcs = get_bcs(V, Time)

            # solve problem
            solve(F == 0, U_kh, bcs)

            solutions_v = [U_kh.sub(i, deepcopy=True).vector() for i in range(Time.n_dofs)]
            offset_v = Time.n_dofs
            solutions_p = [U_kh.sub(i + offset_v, deepcopy=True).vector() for i in range(Time.n_dofs)]
            solutions_e = [U_kh.sub(i + 2*offset_v, deepcopy=True).vector() for i in range(Time.n_dofs)]


            # get v0 for next slab
            U0.vector().set_local(np.concatenate((
                Time.get_solution_at_time(slab[1]-Time.epsilon, solutions_v),
                Time.get_solution_at_time(slab[1]-Time.epsilon, solutions_p), #np.zeros((self.Vh.sub(1).dim(),))
                Time.get_solution_at_time(slab[1]-Time.epsilon, solutions_e)
            )))
            v0, p0, e0 = split(U0)
    
            # plot final solution on slab
            print(f"t = {slab[1]}:")
            vfile << (U0.split(deepcopy=True)[0], slab[1])
            efile << (U0.split(deepcopy=True)[2], slab[1])
            
            if self.goal_functional is not None:
                self.func_vals.append(self.goal_functional(Time, v0, p0, e0))
            print("Done.\n")
        


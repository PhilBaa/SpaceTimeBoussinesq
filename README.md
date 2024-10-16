# Space-Time FEM for Boussinesq equations
Tensor-product space-time Finite Element Method for the Boussinesq equations in FEniCS

Group project for the course "Space-time methods" at the Leibniz University Hannover in the winter semester 2023/2024 <br>
Authors: [Philipp Baasch](https://github.com/PhilBaa), [Göran Ratz](https://github.com/GARHG)  
Tutors: [Julian Roth](https://github.com/mathmerizing/), [Hendrik Fischer](https://github.com/Hendrik240298)

The bulk of the programm is contained in solver.py. schaefer_turek, rectangles and laser implement (modified) problems found in [[1]](#1),[[2]](#2) and [[3]](#3) respectively.

## Strong and Weak Formulation

The Boussinesq equations are essentially a coupled version of the Navier-Stokes and heat equation. For an incompressible fluid with velocity field $v$ and pressure field $p$, temperature field $\theta$ and space-time domain $\Omega = \Omega_s \times I$ the strong form is [[4]](#4)

$$
  \rho\frac{\partial v}{\partial t}-\nabla\cdot (\rho\hspace{0.1cm}\eta(\theta)\hspace{0.1cm}\nabla v)+\nabla p+(\rho v\cdot\nabla)\hspace{0.1cm}v-\alpha\theta g=0\quad \text{in}\quad \Omega,
$$

$$
  \nabla\cdot v=0 \quad \text{in}\quad \Omega,
$$

$$
  \frac{\partial \theta}{\partial t}-\nabla\cdot(k\hspace{0.1cm} \nabla\theta)+v\cdot\nabla\theta=f \quad \text{in}\quad \Omega,
$$

with suitable boundary conditions, which depend on the problem. $\rho$ and $k$ are the density and thermal conductivity respectively, while $\mathcal{f}$ describes a heat source. The last term in the first equation $-\alpha\theta g$ models the buoyancy force due to a gravitational field $g$ and thermal expansion coefficient $\alpha$. In this formulation the viscosity is temperature dependent according to the Arrhenius equation

$$
  \eta(\theta):=\eta_0\exp{\big(\frac{E_A}{R(\theta + T_0)}\big)}.
$$

Assuming Dirchlet data $v_D$ and $\theta_D$, we introduce function spaces

$$
  X_\theta:=\theta_D+\lbrace\theta\in L^2(I,H^1_0(\Omega)): \partial_t \theta \in L^2(I,(H^{1}_{0}(\Omega))^*)\rbrace,
$$

$$
  X_v:=v_D+\lbrace\theta\in L^2(I,H^1_0(\Omega)^d): \partial_t v \in L^{2}(I,(H^1_0(\Omega)^d)^*)\rbrace,
$$

$$
  X_p:= L^{2}(I,L^0(\Omega)),
$$

where $L^0 := \lbrace p \in L^2(\Omega) : \int_\Omega p\, dx = 0\rbrace $. Using Taylor-Hood elements in space and dG(r) in time, we obtain the weak formulation:

Find $u = (v, p, \theta) \in X_v \times X_p \times X_v$ such that

$$
  a(u, \phi) = \sum^M_{m=1}\int_{t_{m-1}}^{t_m}[(\rho\hspace{0.1cm}\partial_t v, \phi^v)+(\partial_t\theta, \phi^\theta)+(\eta(\theta)\hspace{0.1cm}\nabla v,\nabla \phi^{v})+(\rho v\cdot\nabla v,\phi^v)
$$

$$
  -(p,\nabla\phi^v)+(\nabla\cdot v,\phi^p)+(k(\theta)\hspace{0.1cm}\nabla\theta,\nabla\phi^{\theta}))+(v\cdot\nabla\theta, \phi^{\theta})+(\alpha g\theta,\phi^v)\big]\mathrm{dt}
$$

$$
  +\sum^{M-1}_{m=1}(\rho[v]_m,\phi^{v,+}_m)+([\theta]_m,\phi^{\theta,+}_m)+(\rho v^+_0,\phi^{v,+}_0)+(\theta^+_0,\phi^{\theta,+}_0) 
$$

$$
  =\sum^M_{m=1}\int_{t_{m-1}}^{t_m} (f, \phi^\theta)+(v^0, \phi^{v,+}_0)+(\theta^0,\phi^{\theta}_0)=:F^{dG}(\phi)\quad\forall \phi = (\phi^v, \phi^p, \phi^\theta) \in X_v \times X_p \times X_v,
$$

where $(\cdot, \cdot)$ is the usual $L^2$ inner product over space. Moreover, $v^0$ and $\theta^0$ are initial conditions. In practice it is computanionally expensive to use this weak form as it is, which is why we divided the time interval into slabs and performed a time marching scheme using dG(1).

## Results
For the Schäfer Turek benchmark [[1]](#1) with additional boundary conditions $\theta = 600$ at $\partial\Omega_\text{Cylinder}$ and $\theta=0$ at the other boundaries, and parameters $R=E_A$, $T_0$, $k=0.005$, $g = (0, -9.81)$ and $\alpha=0.005$, we obtain the following results as well as drag and lift values:

https://github.com/user-attachments/assets/40fe9dae-a3dd-4a57-9686-6f5a01437ec5

![drag_lift](https://github.com/user-attachments/assets/c188a412-8c4b-47d9-80af-88a831d05de0)


## References


<a id="1">[1]</a> 
Schäfer, M. & Turek, S. & Durst, F. & Krause, E. & Rannacher, R. (1996). Benchmark Computations of Laminar Flow Around a Cylinder. [DOI: 10.1007/978-3-322-89849-4_39](https://doi.org/10.1007/978-3-322-89849-4_39)

<a id="2">[2]</a> 
Cai, S. & Wang, Z. & Wang, S. & Perdikaris, P. & and Karniadakis, G. E. (2021). Physics-Informed Neural Networks for Heat Transfer Problems.  [DOI: 10.1115/1.4050542](https://doi.org/10.1115/1.4050542)

<a id="3">[3]</a> 
Beuchler, S. & Endtmayer, B.  & Lankeit, J. & Wick, T. (2024). Multigoal-oriented a posteriori error control for heated material processing using a generalized Boussinesq model. [DOI: 10.5802/crmeca.160](https://doi.org/10.5802/crmeca.160 )

<a id="4">[4]</a> 
Thiele, Jan Philipp. Doctoral dissertation at Leibniz University Hannover.

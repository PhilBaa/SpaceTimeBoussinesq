# Space-Time FEM for Boussinesq equations
Tensor-product space-time Finite Element Method for the Boussinesq equations in FEniCS

Group project for the course "Space-time methods" at the Leibniz University Hannover in the winter semester 2023/2024 <br>
Authors: [Philipp Baasch](https://github.com/PhilBaa), [Göran Ratz](https://github.com/GARHG)  
Tutors: [Julian Roth](https://github.com/mathmerizing/), [Hendrik Fischer](https://github.com/Hendrik240298)



## Strong and Weak Formulation

The Boussinesq equations are essentially a coupled verion of the Navier-Stokes and heat equation. For an incompressible fluid with velocity field $v$ and pressure field $p$, temperature field $\theta$ and space-time domain $\Omega = \Omega_s \times I$ the strong form is 

$$
  \rho\frac{\partial v}{\partial t}-\nabla\cdot (\rho\nu(\theta)\nabla v)+\nabla p+(\rho v\cdot\nabla)v-\alpha\theta g=0\quad \text{in}\quad \Omega,
$$

$$
  \nabla\cdot v=0 \quad \text{in}\quad \Omega,
$$

$$
  \frac{\partial \theta}{\partial t}-\nabla\cdot(k \nabla\theta)+v\cdot\nabla\theta=f \quad \text{in}\quad \Omega
$$

with suitable boundary conditions, which depend on the problem. $\rho$ and $k$ are the density and thermal conductivity respectively. f describes a heat source. The last term in the first equation $-\alpha\theta g$ models the buoyancy force due to a gravitational field $g$ and thermal expansion coefficient $\alpha$. In this formulation the viscosity is temperature dependent according to the Arrhenius equation

$$
  \nu(\theta):=\nu_0\exp{\frac{E_A}{R(\theta + T_0)}}.
$$

Assuming Dirchlet data $v_D$ and $\theta_D$ we introduce function spaces

$$
  X_\theta:=\theta_D+\{\theta\in L^2(I,H^1_0(\Omega)): \partial_t \theta \in L^2(I,(H^{1}_{0}(\Omega))^*\},
$$

$$
  X_v:=v_D+\{\theta\in L^2(I,H^1_0(\Omega)^d): \partial_t v \in L^{2}(I,(H^1_0(\Omega)^d)^*)\},
$$

$$
  X_p:= L^{2}(I,L^0(\Omega)),
$$

where $L^0 := \{p \in L^2(\Omega) : \int_\Omega p\, dx = 0\}$. Using Taylor-Hood elements and space and dG(r) in time we obtain the weak formulation: Find $u = (v, p, \theta) \in X_v \times X_p \times X_v$ such that

$$
  a(u, \phi) = \sum^M_{m=1}\int_{t_{m-1}}^{t_m}[(\rho\partial_t v, \phi^v)+(\partial_t\theta, \phi^\theta)+(\nu(\theta)\nabla v,\nabla \phi^{v})+(\rho v\cdot\nabla v,\phi^v)
$$

$$
  -(p,\nabla\phi^v)+(\nabla\cdot v,\phi^p)+(k(\theta)\nabla\theta,\nabla\phi^{\theta}))+(v\cdot\nabla\theta, \phi^{\theta})+(\alpha g\theta,\phi^v)\big]\mathrm{dt}
$$

$$
  +\sum^{M-1}_{m=1}(\rho[v]_m,\phi^{v,+}_m)+([\theta]_m,\phi^{\theta,+}_m)+(\rho v^+_0,\phi^{v,+}_0)+(\theta^+_0,\phi^{\theta,+}_0) 
$$

$$
  =\sum^M_{m=1}\int_{t_{m-1}}^{t_m} (f, \phi^\theta)+(v^0, \phi^{v,+}_0)+(\theta^0,\phi^{\theta}_0)=:F^{dG}(\phi)\quad\forall \phi = (\phi^v, \phi^p, \phi^\theta) \in X_v \times X_p \times X_v,
$$

where $(\cdot, \cdot)$ is the usual $L^2$ inner product over space. $v^0$ and $\theta^0$ are initial conditions. In Practice it is computanionally expensive to use this weak form as is, which is why we divided the time interval into slabs and performed a time marching scheme using dG(1).

## Results



## References


<a id="1">[1]</a> 
Beuchler, Sven & Endtmayer, Bernhard  & Lankeit, Johannes & Wick, Thomas. (2024). Multigoal-oriented a posteriori error control for heated material processing using a generalized Boussinesq model. [10.5802/crmeca.160](https://doi.org/10.5802/crmeca.160 )

<a id="2">[2]</a> 
Thiele, Jan Philipp. Doctoral dissertation at Leibniz University Hannover.
Shengze Cai, Zhicheng Wang, Sifan 

<a id="3">[3]</a> 
Cai, Shengze& Wang, Zhicheng & Wang, Sifan & Perdikaris, Paris & and Karniadakis, George Em (2021). Physics-Informed Neural Networks for Heat Transfer Problems.  [10.1115/1.4050542](https://doi.org/10.1115/1.4050542)

<a id="4">[4]</a> 

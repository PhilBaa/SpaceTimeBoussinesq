# Space-Time FEM for Boussinesq equations
Tensor-product space-time Finite Element Method for the Boussinesq equations in FEniCS

Group project for the course "Space-time methods" at the Leibniz University Hannover in the winter semester 2023/2024 <br>
Authors: [Philipp Baasch](https://github.com/PhilBaa), Göran Ratz  
Tutors: [Julian Roth](https://github.com/mathmerizing/), [Hendrik Fischer](https://github.com/Hendrik240298)

## ToDo list
- [ ] Read Sections 2.1 to 2.3 and the problem description in Section 4.1 in [[Roth et al. 2023]](https://doi.org/10.1515/cmam-2022-0200).
- [ ] Look again carefully at the [Navier-Stokes code from Exercise 4](https://github.com/mathmerizing/SpaceTimeFEM_2023-2024/blob/main/Exercise4/Exercise_4_NavierStokes.ipynb), try to understand it and write down your open questions. 
- [ ] Try to play around with the Navier-Stokes code a little:
  - [ ] Compare dG(0) and dG(1) time discretizations. Can you imagine why dG(1) might be more suitable for Navier-Stokes?
  - [ ] Compare $s_v = 1, s_p = 1$ results with $s_v = 2, s_p = 1$. Do you make any strange observations?

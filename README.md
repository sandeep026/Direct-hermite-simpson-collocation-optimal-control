# **Hermite simpson collocation-optimal-control**

An optimal control problem for a cart system is solved using Hermite simpson collocation. The optimal control problem is discretized and transformed to a nonlinear optimization problem and solved using the state of the art solver IPOPT.

# Requirements
- MATLAB/[OCTAVE](https://octave.org/)
- [Casadi](https://web.casadi.org/)

# Problem description

The optimal control problem for the cart system [^1] is provided below. z <sub>1</sub> and z <sub>2</sub> are the position and velocity of the cart and they comprise the states. f is the force applied and there is a drag force which is proportional to the velocity of the cart. The system starts from rest and additional boundary condition is placed at the end of the trajectory. Along the trajectory the control effort is minimized from time 0 to 2.

![image](https://user-images.githubusercontent.com/16457676/236567436-9d87b891-e74f-4299-802c-a394693c1f60.png)

# Analytical solution

The system admits the following analytical solution, which can be later used to verify the numerical solution and its accuracy.

![image](https://user-images.githubusercontent.com/16457676/236629178-b6da4837-b1d8-454d-9ec4-2d67fb1abeba.png)

# Discretization

The time domain is discretized and a medium order polynomial is used to approximate the state and control. The integrals are approximated using the Simpson quadrature and the nonlinear constraints are only evaluated at the grid points. Collectively, this results in a nonlinear optimization problem which can be solved using a off-the-shelf solver. Here, we have used IPOPT for finding the optimal solution. 

# Results

The results are plotted in phase space for a grid size of 50 and they are in close agreement with the analytical solution.

![image](https://github.com/sandeep026/hermite-simpson-collocation--optimal-control/assets/16457676/4ee32144-1d7c-42b2-8d33-6d8dcf4fe342)

![image](https://github.com/sandeep026/hermite-simpson-collocation--optimal-control/assets/16457676/7dfb1faf-06d2-43aa-a7a9-2015f52d8899)

# Known issues


# References

[^1]: Conway, B. A. and K. Larson (1998). Collocation versus differential inclusion in direct optimization. Journal of Guidance, Control, and Dynamics, 21(5), 780–785

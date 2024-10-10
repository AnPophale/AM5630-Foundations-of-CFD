# AM5630-Foundations-of-CFD

This course included two assignments, implementing the Finite Volume Method in MATLAB: 2D diffusion and 2D convection-diffusion for scalar transport. The codes can be found in the [MATLAB codes](https://github.com/AnPophale/AM5630-Foundations-of-CFD/tree/main/MATLAB%20Codes) folder, along with the [reports](https://github.com/AnPophale/AM5630-Foundations-of-CFD/tree/main/Reports). Below is an explanation of the governing equations, discretization schemes, and results for both cases.

### 2D Diffusion:  
We consider the problem of 2D heat conduction as given below,
<p align="center">
  <img src="https://github.com/user-attachments/assets/c73f4ae1-7728-4de1-8926-f8fa3cbec690" alt=Figure 1: Domain for 2D heat conduction" style="width: 50%;">
</p>
<p align="center">
  <em>Figure 1: Domain for 2D heat conduction</em>
</p>

The boundary conditions are,
- **Boundary 1**: $`T = 15 \, ^\circ \text{C}`$
- **Boundary 2**: $`T(y) = 5\left(1 - \frac{y}{H}\right) + 15 \sin\left(\frac{\pi y}{H}\right)`$
- **Boundary 3**: $`T = 10 \, ^\circ \text{C}`$
- **Boundary 4** (Insulated): $`\frac{\partial T}{\partial x} = 0`$

The thermal conductivity is a function of space given as $`k = 16(\frac{y}{H} + 1) `$ and a heat source term is given as $`S = -1.5`$ per unit area.

The governing equation for 2D heat conduction with a source term is given as,
```math
\frac{\partial}{\partial x} \left( k \frac{\partial T}{\partial x} \right) 
+ \frac{\partial}{\partial y} \left( k \frac{\partial T}{\partial y} \right) + S = 0
```

We use the Finite Volume Method with central difference scheme and source term linearization to discretize the equations. In finite volume method, the governing differential equation is first integrated over a control volume. The complete details of the theory for the Finite Volume Method can be found in [1]. Consider a control volume shown in Fig. 2.

<p align="center">
  <img src="https://github.com/user-attachments/assets/78aef86c-a432-41ee-8134-2b53732948c8"  alt=Figure 2: 2D Control Volume for the Finite Volume Method" style="width: 30%;">
</p>
<p align="center">
  <em>Figure 2: 2D Control Volume for the Finite Volume Method</em>
</p>

Integrating the governing differential equation over this 2D control volume, we get,

```math
0 = \int_w^e \int_s^n \frac{\partial}{\partial x} \left( k \frac{\partial T}{\partial x} \right) \, dx \, dy 
+ \int_s^n \int_w^e \frac{\partial}{\partial y} \left( k \frac{\partial T}{\partial y} \right) \, dx \, dy + \int_w^e \int_s^n S \, dx \, dy
```

We assume that the source term remains constant over the volume and hence it is not integrated. Integrating the other terms, 
```math
\left[ k_{e} \Delta y \left( \frac{\partial \phi}{\partial x} \right)_{e} - \, k_{w} \Delta y \left( \frac{\partial \phi}{\partial x} \right)_{w} \right] + \left[ k_{n} \Delta x \left( \frac{\partial \phi}{\partial y} \right)_{n} - \, k_{s} \Delta x \left( \frac{\partial \phi}{\partial y} \right)_{s} \right] + S \Delta x \Delta y = 0
```

All the variables are only stored at the cell centers, hence to calculate the gradients at the cell faces, we need different schemes. The central difference scheme uses following equation to calculate the fluxes at the faces using the values at the cell centers.

```math
\left( \frac{dT}{dx} \right)_{e} = \frac{T_{E} - T_{P}}{\delta x_e}, \quad 
\left( \frac{dT}{dx} \right)_{w} = \frac{T_{P} - T_{W}}{\delta x_w}, \quad 
\left( \frac{dT}{dx} \right)_{n} = \frac{T_{N} - T_{P}}{\delta y_n}, \quad 
\left( \frac{dT}{dx} \right)_{s} = \frac{T_{P} - T_{S}}{\delta y_s}
```

Substituting this into the previous equation and rearranging it in a form given as $`a_P T_P = a_E T_E + a_W T_W + a_N T_N + a_S T_S + S \Delta x \Delta y `$ we have
```math
a_{W} = \frac{k_{w} \Delta y}{\delta x_{e}}, \quad 
a_{E} = \frac{k_{e} \Delta y}{\delta x_{w}}, \quad 
a_{S} = \frac{k_{s} \Delta x}{\delta x_{n}}, \quad 
a_{N} = \frac{k_{n} \Delta x}{\delta x_{s}}, \quad 
a_{P} = a_{W} + a_{E} + a_{S} + a_{N} 
```

A negative constant source term can lead to diverging results during the numerical procedure to solve the system of equations, hence we make the following adjustment for numerical stability shifting the source term to the LHS from the RHS, using $`T_old`$ which is the value of the temperature from the previous iteration.
```math
a_{P} = a_{W} + a_{E} + a_{S} + a_{N} + \frac{S \Delta x \Delta y}{T_{old}}
```

These equations are applied to each control volume in the discretized domain and the resulting system of linear equations is solved using the Gauss Seidel method. The code implementing the above described method is given [here](https://github.com/AnPophale/AM5630-Foundations-of-CFD/blob/main/MATLAB%20Codes/2D_Diffusion.m).

**Results**:  
The temperature distribution in the domain using a 80x80 mesh is shown below

<p align="center">
  <img src="https://github.com/user-attachments/assets/0e842e6d-f90d-4c3d-bcb2-944f632b78a7"  alt="Figure 3: Temperature field in the domain" style="width: 50%;">
</p>
<p align="center">
  <em>Figure 3: Temperature field in the domain</em>
</p>

A grid convergence study is also performed to find the optimal mesh size. A stretch mesh with smaller elements near the boundary with large gradients is also used and the results are compared with larger grid sizes to compare the performance. The variation of the residual vs number of iteration for various convergence criteria is also conducted. The effect of changing the boundary conditions is also explored and the details of all the analysis can be found in the [report](https://github.com/AnPophale/AM5630-Foundations-of-CFD/blob/main/Reports/2D%20Diffusion%20Report.pdf).

### 2D Convection-Diffusion:  
Here, we consider the problem of heat convection and conduction with a given velocity field as shown in Fig. 4 below along with given mesh data.
<p align="center">
  <img src="https://github.com/user-attachments/assets/56bc87d1-6d62-448b-829b-838403af134e" alt="Figure 4: Computational domain and the given velocity field" style="width: 80%;">
</p>
<p align="center">
  <em>Figure 4: Computational domain and the given velocity field</em>
</p>

The dimensions for the domain and given physical values are $\rho = 1$, $k/C_p = 1/50$, and $h_A/H = h_C/H = 0.068$. The boundary conditions are $U_A = 1$, $U_B = 0$, $U_C = 1$, $V_D = 0$, $T_A = 20^\circ C$, and at $x = L$ (other than the outlet), $T = 50^\circ C$. The governing equation for the temperature with convection and diffusion is given as,
```math
\frac{\partial}{\partial x} (\rho U T) + \frac{\partial}{\partial y} (\rho V T) = \frac{\partial}{\partial x} \left( \Gamma \frac{\partial T}{\partial x} \right) + \frac{\partial}{\partial y} \left( \Gamma \frac{\partial T}{\partial y} \right) + S, \quad \text{where } \Gamma = \frac{k}{C_p}
```

We apply the FVM discretization procedure as described previously using the central difference scheme for diffusive terms but here, we use the hybrid scheme for the convective term which is explained below.
On integration, the convective term gives
```math
(\rho U T)_{w} \Delta y - (\rho U T)_{e} \Delta y + (\rho V T)_{n} \Delta x - (\rho V T)_{s} \Delta x
```
Again as done in the FVM, all variables are stored at cell centers and values are not known at the cell faces. As the velocity profile is already known at the cell centers, it can be interpolated to the cell faces. But the calculate the unknown T at the cell faces, we need some discretization scheme for which we use the hybrid scheme. 

The temperature values at the cell faces is calculated in the hybrid scheme as an upwind or central differencing scheme based on the cell Peclet number.  Essentially, based on the ratio of convection to diffusion at each cell, either the central difference scheme (suitable for diffusion) or the first order upwind scheme (suitable for convection) is used. The formulation for the hybrid scheme is given considering the first term ($`T_w`$) in the above equation and it is written as,

```math
T_w = 
\begin{cases} 
\frac{T_W + T_P}{2} & \text{if } |Pe_w| < 2 \quad (\text{Central Difference}) \\ 
T_W & \text{if } Pe_w > 2 \quad (\text{First order Upwind}) \\ 
T_P & \text{if } Pe_w < -2 \quad (\text{First order Upwind}) 
\end{cases}
\quad \text{where } Pe_w = \frac{\rho U_w \Delta x_w}{\Gamma}
```

This is applied to each face to get the discretized equations. The system of linear equations generated is solved using the Gauss Seidel as well as an iterative 2D Tri Diagonal Matrix Algorithm (TDMA) considering different sweep directions. The MATLAB code for this is given in the codes folder as 2D_Convection_Diffusion.m 

**Results:** 

**References:**  
[1] H. Versteeg and W. Malalasekera. An Introduction to Computational Fluid Dynamics - The Finite Volume Method. Longman Scientific & Technical, Harlow, England, 1st edition, 1995.

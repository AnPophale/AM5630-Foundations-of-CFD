# AM5630-Foundations-of-CFD

This course included 2 assignments using the Finite Volume Method - 2D diffusion and 2D convection-diffusion for scalar transport. The codes can be found in the MATLAB codes folder along with the reports. Following is an explanation of the governing equations, discretization schemes and results for both the cases.

### 2D Diffusion:  
We consider the problem of 2D heat conduction as given below

![image](https://github.com/user-attachments/assets/c73f4ae1-7728-4de1-8926-f8fa3cbec690)


<p align="center">
  <img src="https://github.com/user-attachments/assets/959fec8f-649d-48de-ac8b-84a3b29b223c" alt="Comparison of turbulence kinetic energy using Reynolds Stress model with wall function and DNS data" style="width: 50%;">
</p>
<p align="center">
  <em>Figure 4: Comparison of turbulence kinetic energy using Reynolds Stress model with wall function and DNS data</em>
</p>

The boundary conditions are 

The thermal condictivity is a function of space given as
The governinng equation is given as 

Where the source term is given as 

We use the Finite Volume Method with central difference scheme and source term linearization to discretize the equations. In finite volume method, the governing differential equation is first integrated over a control volume. Consider a control volume given as below. Integrating the equation over it, we get

We assume that the source term remains constant over the volume and hence it is not integrated

All the variables are only stored at the cell centers, hence to calculate the values at the cell faces, we need different schemes. The central differnce scheme uses following equation to calculte the fluxes at the faces using the values at the cell centers

Substituting this and simplifying, we get the follwing equation

A negative constant source term can give divergence during the numerical procdedure to solve the system of equations, hence we make the follwing adjustment for numerical stability

These equations are applied to each control volume in the discretized domain and the resulting system of linear equations is solved using the Gauss Seidel method

**Results**:  

### 2D Convection-Diffusion:  




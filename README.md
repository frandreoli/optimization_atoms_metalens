# Introduction

In \[[1](#Andreoli2023b)\], the concept was introduced of a thin metalens composed only of arrays of atoms, with spatially varying lattice constants. The design of such an atomic metalens depends on three free parameters, that are the disk thickness $0<\Delta R<1$, the overall phase shift $-\pi<\phi\_0\leq \pi$ and the buffer fraction $0\leq \alpha\leq 0.5$. Additional details of the physical meaning of these variables are provided at Ref. \[[1](#Andreoli2023b)\]. 

The metalens efficiency $\eta$ can be numerically computed via exact coupled-dipole equations, as further implemented and described at the repository [frandreoli/atoms_optical_response](https://github.com/frandreoli/atoms_optical_response). Maximizing $\eta$ as a function of the free physical variable is a global optimization problem with bounded parameters. The remarkably small number of variables and their rectangular bounds highly reduce the parameter space, largely helping to simplify the optimization. On the other hand, the object function is a highly irregular (and in principle non-continuous and non-differentiable) function, which exhibits a large number of local minima. 



# Options

## Solvers

# References 

<a id="Andreoli2023b">[1]</a> 
Andreoli F, High AA, Chang DE, 
*Metalens formed by structured, sub-wavelength atomic arrays*, 
paper in preparation (2023)

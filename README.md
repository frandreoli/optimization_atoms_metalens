# Introduction

In Ref. \[[1](#Andreoli2023b)\], the concept was introduced of a thin metalens composed only of arrays of atoms, with spatially varying lattice constants. The design of such an atomic metalens depends on three free parameters, that are the disk thickness $0<\Delta R<1$, the overall phase shift $-\pi<\phi\_0\leq \pi$ and the buffer fraction $0\leq \alpha\leq 0.5$. Additional details of the physical meaning of these variables are provided at Ref. \[[1](#Andreoli2023b)\] and in the description of the repository  [_frandreoli/atoms_optical_response_](https://github.com/frandreoli/atoms_optical_response). 

## Overview of the problem

The metalens efficiency $\eta$ can be numerically computed via exact coupled-dipole equations, as further implemented and described in the repository [_frandreoli/atoms_optical_response_](https://github.com/frandreoli/atoms_optical_response). The computational complexity of this core operation roughly scales as the creation and inversion of a $N\times N$, complex-symmetric matrix, where $N$ is the number of atoms. Maximizing $\eta(\Delta R, \phi\_0, \alpha)$ is a global optimization problem over three bounded parameters. The remarkably small number of variables and their rectangular bounds highly reduce the parameter space, largely helping to simplify the optimization. 

On the other hand, the object function is highly irregular (and in principle non-continuous and non-differentiable), and it exhibits an extremely large number of local minima. This is further complicated by the dependence of $N$ from the set of variables $\Delta R$, $\phi\_0$ and $\alpha$. Overall, relevant sizes of the metalens are associated to number of atoms as large as $N\sim 10^4-10^5$, which can make the evaluation of $\eta(\Delta R, \phi\_0, \alpha)$ remarkably slow, even when speeding up the linear algebra by feeding (up to) $32$ _multiple cores_ to OpenBlas or accelerating other internal operations by parrallelizing them on $\geq 32$ cores. 



# Options
The definitions of the main physical variables and some options (namely, `z_fixed_option` and `phase_center_ring_option`) is defined in the repository [_frandreoli/atoms_optical_response_](https://github.com/frandreoli/atoms_optical_response). Additional options are available for this specific optimization problem, that we describe below.

### System options
- `fill_until_r_lens_option` \
When set to `true`, atoms are positioned up to the fixed radius value `r_lens`, even if this means the introduction of a fraction of the last ring. Otherwise, when set to `false`, then only integer rings are accepted, meaning that the effective radius of the lens can be $\leq$`r_lens`. This latter choice makes the number of atoms $N$ and efficiency $\eta$ depend more sharply on the optimization variables, and it is discouraged.

## Solvers

# References 

<a id="Andreoli2023b">[1]</a> 
Andreoli F, High AA, Chang DE, 
*Metalens formed by structured, sub-wavelength atomic arrays*, 
paper in preparation (2023)

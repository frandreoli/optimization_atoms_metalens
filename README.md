# Introduction

In Ref. \[[1](#Andreoli2023b)\], the concept was introduced of a thin metalens composed only of arrays of atoms, with spatially varying lattice constants. The design of such an atomic metalens depends on three free parameters, that are the disk thickness $0<\Delta R<1$, the overall phase shift $-\pi<\phi\_0\leq \pi$ and the buffer fraction $0\leq \alpha\leq 0.5$. Additional details of the physical meaning of these variables are provided at Ref. \[[1](#Andreoli2023b)\] and in the description of the repository  [_frandreoli/atoms_optical_response_](https://github.com/frandreoli/atoms_optical_response). 

## Overview of the problem

The metalens efficiency $\eta$ can be numerically computed via exact coupled-dipole equations, as further implemented and described in the repository [_frandreoli/atoms_optical_response_](https://github.com/frandreoli/atoms_optical_response). The computational complexity of this core operation roughly scales as the creation and inversion of a $N\times N$, complex-symmetric matrix, where $N$ is the number of atoms. Maximizing $\eta(\Delta R, \phi\_0, \alpha)$ is a global optimization problem over three box-bounded parameters. The remarkably small number of variables and their rectangular bounds highly reduce the parameter space, largely helping to simplify the optimization. 

On the other hand, the object function is highly irregular (and in principle non-continuous and non-differentiable), and it exhibits an extremely large number of local minima. This is further complicated by the dependence of $N$ from the set of variables $\Delta R$, $\phi\_0$ and $\alpha$. Overall, relevant sizes of the metalens are associated to number of atoms as large as $N\sim 10^4-10^5$, which can make the evaluation of $\eta(\Delta R, \phi\_0, \alpha)$ remarkably slow, even when speeding up the linear algebra by feeding (up to) $32$ _multiple cores_ to OpenBlas or accelerating other internal operations by parrallelizing them on $\geq 32$ cores. 



# Setting the code
The definitions of the main physical variables and some options (namely, `z_fixed_option` and `phase_center_ring_option`) is defined in the repository [_frandreoli/atoms_optical_response_](https://github.com/frandreoli/atoms_optical_response). Additional options are available for this specific optimization problem, that we describe below. The result of the optimization are printed at each step into the `stdout` (by flushing it at each step, to avoid loss of data after a potential crash). In a linux environment, a _shell_ wrapper is provided to redirect the `stdout` and `stderr` on two different files with proper labels.

## Solvers
The solvers can be chosen by defining a proper value of the integer variable `solver_algorithm_index`, chosen among the following list:

- From [_BlackBoxOptim_](https://github.com/robertfeldt/BlackBoxOptim.jl):
1)  `:adaptive_de_rand_1_bin` _Differential Evolution optimizer (metaheuristics)._
2)  `:adaptive_de_rand_1_bin_radiuslimited` _Differential Evolution optimizer (metaheuristics) with limited radius._ 
3)  `:resampling_memetic_search` _Memetic algorithm._
4)  `:resampling_inheritance_memetic_search` _Memetic algorithm with inheritance._
5)  `:simultaneous_perturbation_stochastic_approximation` _Stochastic Approximation algorithm._
6)  `:separable_nes` _Natural Evolution Strategies - separable._
7)  `:xnes` _Natural Evolution Strategies - exponential._
8)  `:dxnes` _Natural Evolution Strategies - exponential - distance-weighted._
9)  `:generating_set_search` _Generating-set direct search._
10) `:probabilistic_descent` _Generating-set direct search, with probabilistic descent._

- From Optim (https://julianlsolvers.github.io/Optim.jl/stable/):
11) `ParticleSwarm()` _Particle Swarm Optimization_
12) `SAMIN()` _Simulated Annealing with intrinsic bounds_
13) `SimulatedAnnealing()` _Simulated Annealing*_
14) `NelderMead()` _Nelder-Mead*_

_*All the solvers are design to directly deal with box-bounded problems, except for choices 13 and 14, where the box constraint is embedded inside the object function (i.e. this latter is set to return a very inefficient value when the parameters exit the box bounds. This is generally not recommendable, and their behaviour is indeed underwhelming._ 

## Optimization settings
Specific settings related to the chosen solver can be implemented when calling the solver in the code. Generic settings are the following:

- `thickness_range, phase_range, buffer_range` \
Choice of the box ranges of accepted parameters, that are automatically set to: thickness $0.1<\Delta R<0.8$, phase $-\pi<\phi\_0\leq \pi$ and buffer $0\leq \alpha\leq 0.5$.

- `max_steps, max_time_sec, max_stuck_sequence` \
Maximum allowed amount of, respectively, optimization steps, computational time (in seconds) and consecutive steps where the algorithm is stuck with the same solution.

## Options
Global options for the code are available, which are described here below.
### System options
- `fill_until_r_lens_option` \
If set to `true`, then atoms are positioned up to the fixed radius value `r_lens`, even if this means the introduction of a fraction of the last ring. Otherwise, when set to `false`, then only integer rings are accepted, meaning that the effective radius of the lens can be $\leq$`r_lens`. This latter choice makes the number of atoms $N$ and efficiency $\eta$ depend more sharply on the optimization variables, and it is discouraged.

### Optimization options
- `initial_guess_option` \
If set to `true`, then an initial guess `initial_guess` is fed into the chosen optimization algorithm. Otherwise, when set to `false`, then `initial_guess` is randomly initialized.

- `monotonic_escape_option` \
If set to `true`, then the code stops the optimization flow if it detects that a worst efficiency is proposed as the new optimal one. We observed this unexpected behaviour with the solver `SAMIN()` of the library [_Optim_](https://github.com/JuliaNLSolvers/Optim.jl).

### Supplementary options
- `debug_r_atoms_option` \
If set to `true`, then the code does not perform any optimization, but it rather export on an `HDF5` files the exact atomic positions of an illustrative atomic metalens.


# References 

<a id="Andreoli2023b">[1]</a> 
Andreoli F, High AA, Chang DE, 
*Metalens formed by structured, sub-wavelength atomic arrays*, 
paper in preparation (2023)

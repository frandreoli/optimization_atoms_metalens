#Designed for Julia > 1.5
using LinearAlgebra, Dates, HDF5, Random, Distributed
include("Metalens OPTIM - SM Functions.jl")
include("Metalens OPTIM - Core Evaluation.jl")
#


################## OPTIONS ####################################################################################################################################
#
const z_fixed_option               =      ["YES" ; "NO"][1] #must be yes
const phase_center_ring_option     =      ["YES" ; "NO"][1] #must be yes
const initial_guess_option         =      [true ; false][1]
const use_blackboxoptim_option     =      [true ; false][2]
#


################## FIXED PARAMETERS ####################################################################################################################################
#
Random.seed!()
if nworkers()==1 addprocs() end
BLAS.set_num_threads(nworkers())
const TN = Float32
const nBulk = 1.0
const lambda0 = 1.0/nBulk
const k0 = 2.0*pi/lambda0
const dipoles_polarization = [1.0 ; 0.0 ; 0.0]
const field_polarization   = [1.0 ; 0.0 ; 0.0]
length(ARGS)>=1 ? args_checked=ARGS[:] : args_checked=["_vacuum","_mirror"][1:1]
#


################## GEOMETRICAL PARAMETERS ##############################################################################################################################
#
const w0                =    4.0*lambda0
const focal_point       =    20*lambda0
const r_lens            =    6*lambda0
const gamma_prime       =    5.75
const laser_detuning    =    0.0
const w0_x              =    w0
const w0_y              =    w0
#
println("")
println("#"^25)
println("")
println("The initial settings are: ")
println("w0 = ",w0)
println("gamma_prime = ", gamma_prime)
println("r_lens = ", r_lens)
println("focal_point = ", focal_point)
println("")
println("#"^25)
#


################## OPTIMIZATION ###########################################################################################################
#
#Optimization parameters
thickness_range = (0.1,0.8)
phase_range = (-3.141592653589793,3.141592653589793)
buffer_range = (0.,0.5)
initial_guess = [2.0/3, 1.0775,0.20048828125]#0.2]
#
#Definition of the solver
if use_blackboxoptim_option
    #
    using BlackBoxOptim
    chosen_solver = (
        :adaptive_de_rand_1_bin,                             #1  Differential Evolution optimizer (metaheuristics)
        :adaptive_de_rand_1_bin_radiuslimited,               #2  Suggested by developers
        :resampling_memetic_search,                          #3  Memetic algorithm 
        :resampling_inheritance_memetic_search,              #4  Memetic algorithm #VERY BAD
        :simultaneous_perturbation_stochastic_approximation, #5  Stochastic Approximation algorithm #ERROR
        :separable_nes,                                      #6  Natural Evolution Strategies
        :xnes,                                               #7  Natural Evolution Strategies  #####
        :dxnes,                                              #8  Natural Evolution Strategies  
        :generating_set_search,                              #9  Generating set (direct) search #BEST UP TO NOW #####
        :probabilistic_descent,                              #10 Generating set (direct) search #####
        :borg_moea,                                          #11 http://borgmoea.org/ #Multi-objective optimization
    )[7]
else
    using Optim
    chosen_solver = (
        ParticleSwarm(; n_particles = 3),                    #1  Particle Swarm Optimization
        SimulatedAnnealing(),                                #2  Simulated Annealing
    )[1]
end
#
#Printing data
println("")
println("The parameter ranges are: ")
println("thickness_range = ", thickness_range)
println("phase_range = ", phase_range)
println("buffer_range = ", buffer_range)
if initial_guess_option
    println("The initial guess is: ")
    println("- thickness guess: ", initial_guess[1])
    println("- phase guess: "    , initial_guess[2])
    println("- buffer guess: "   , initial_guess[3])
end
println("Using the solver: "   , String(Symbol(chosen_solver)))
println("")
println("#"^25)
println("\n\n")
flush(stdout)
#
#Objective function
function objective_func(x) 
    #thickness, phase, buffer
    return 1.0-SM_main(x[1], x[2], x[3])
end
#
#Initializing the optimization
if initial_guess_option
    println("** STARTING the optimization with an initial guess **")
    blackboxoptim_input = (objective_func, initial_guess)
    initial_guess_Optim = initial_guess
else
    println("** STARTING the optimization with NO initial guess **")
    blackboxoptim_input = (objective_func)
    initial_guess_Optim = rand_range.([thickness_range, phase_range, buffer_range])
end
#
function g_func!(G,x)
    error("You shouldn't call the gradient function.")
end
#Optimization core
if use_blackboxoptim_option
    #
    optim_results = bboptimize(blackboxoptim_input...; 
        Method = chosen_solver, 
        SearchRange = [thickness_range, phase_range, buffer_range], 
        NumDimensions = 3,
        TraceMode = :compact, #:verbose, #:silent,
        TraceInterval = 1,#1,
        PopulationSize = 50, #default is 50
        CallbackFunction = x -> begin println("Efficiency = $(1-best_fitness(x)), $(best_candidate(x))"); flush(stdout) end,
        CallbackInterval = 10,#1
        TargetFitness = 0.0
    )
    #
    x_results = best_candidate(optim_results)
    #
else
    #
    lower_range = (x->sort(collect(x))[1]).([thickness_range, phase_range, buffer_range])
    upper_range = (x->sort(collect(x))[2]).([thickness_range, phase_range, buffer_range])
    #
    optim_results = Optim.optimize(objective_func, g_func!, 
        lower_range, 
        upper_range, 
        initial_guess_Optim,
        Fminbox(chosen_solver), 
        Optim.Options(
            show_trace = true, 
            show_every=1, 
            callback = x -> begin flush(stdout) end,
            extended_trace=true
            )
    )
    #
    x_results = optim_results.minimizer
    #
end
#
#Printing the final results
println("#"^25)
println("\n")        
println("** Optimization result: " , best_fitness(optim_results))
println("** Optimal parameters: ")
println("   1) disks_thickness     = "        , x_results[1])
println("   2) phase_shift         = "        , x_results[2])
println("   3) buffer_smooth       = "        , x_results[3])
#


################## FUNCTIONS #############################################################################################################
#
function rand_range(range)
    sorted_range = sort(collect(range))
    (sorted_range[2]-sorted_range[1])*rand()+sorted_range[1]
end

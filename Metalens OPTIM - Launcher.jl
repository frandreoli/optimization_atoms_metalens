#Designed for Julia > 1.5
using LinearAlgebra, Dates, HDF5, Random, Distributed
include("Metalens OPTIM - SM Functions.jl")
include("Metalens OPTIM - Core Evaluation.jl")
#


################## OPTIONS ####################################################################################################################################
#
const z_fixed_option               =      ["YES" ; "NO"][1] #must be yes
const phase_center_ring_option     =      ["YES" ; "NO"][1] #must be yes
#const shift_to_match_option       =      ["YES" ; "NO"][2] #must be no
#const phase_ring_integral_option  =      ["YES" ; "NO"][2] #must be no (?)
const initial_guess_option         =      [true ; false][1]
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
#Defines the filename
#=file_name="r"*string(r_lens/lambda0)[1:min(length(string(r_lens/lambda0)),3)]
file_name*="_f"*string(focal_point/lambda0)[1:min(length(string(focal_point/lambda0)),3)]
file_name*="_w"*string(w0/lambda0)[1:min(length(string(w0_x/lambda0)),3)]
file_name*="_n"*string(nBulk)[1:min(3,length(string(nBulk)))]
file_name*="_gPr"*string(gamma_prime)[1:min(4,length(string(gamma_prime)))]
file_name*="_"*string(Dates.today())
file_name="Metalens.OPTIM_"*file_name*args_checked[1]
#
#Initializes the folder
println("\nTime: ",now(),"\nStarting evaluation of ", @__FILE__,"\n")
println("Output file name: ",file_name,"\n\n")
final_path_name="Data Out/"*file_name*"/"
mkpath(final_path_name)
#
#Writing settings data file
h5write_multiple(final_path_name*"inputs",        [("inputs",        [lambda0 w0 r_lens focal_point nBulk gamma_prime] )])
#
=#
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
#
#Optimization set up
#
using BlackBoxOptim
#
thickness_range = (0.1,0.8)
phase_range = (-pi,pi)
buffer_range = (0,0.5)
initial_guess = [2.0/3, 1.0775,0.20048828125]#0.2]
#
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
#
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
println("Using the solver: "   , chosen_solver)
println("")
println("#"^25)
println("\n\n")
flush(stdout)
#
function objective_func(x)
    # (disks_thickness, phase_shift, buffer_smooth) = Tuple(x)
    return 1.0-SM_main(x[1], x[2], x[3])
end
#
if initial_guess_option
    #
    println("** STARTING the optimization with an initial guess **")
    #
    optim_results = bboptimize(objective_func, initial_guess; 
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
else
    #
    println("** STARTING the optimization with NO initial guess **")
    #
    optim_results = bboptimize(objective_func; 
        Method = chosen_solver, 
        SearchRange = [thickness_range, phase_range, buffer_range], 
        NumDimensions = 3,
        TraceMode = :compact, #:verbose, #:silent,
        TraceInterval = 1,
        PopulationSize = 50, #default is 50
        CallbackFunction = x -> begin println("Efficiency = $(1-best_fitness(x)), $(best_candidate(x))"); flush(stdout) end,
        CallbackInterval = 10,#1
        TargetFitness = 0.0
    )
    #
end
#
x_results=best_candidate(optim_results)
#
println("#"^25)
println("\n")        
println("** Optimization result: " , best_fitness(optim_results))
println("** Optimal parameters: ")
println("   1) disks_thickness     = "        , x_results[1])
println("   2) phase_shift         = "        , x_results[2])
println("   3) buffer_smooth       = "        , x_results[3])
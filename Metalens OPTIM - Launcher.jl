#Designed for Julia > 1.5
using LinearAlgebra, Dates, HDF5, Random, Distributed
using Suppressor
include("Metalens OPTIM - SM Functions.jl")
include("Metalens OPTIM - Core Evaluation.jl")
const intInf = typemax(Int)
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
const r_lens            =    3*lambda0
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
const thickness_range = (0.1,0.8)
const phase_range = (-3.141592653589793,3.141592653589793)
const buffer_range = (0.,0.5)
const initial_guess = [2.0/3, 1.0775,0.20048828125]#0.2]
const max_steps    = [intInf ; 2][2]
const max_time_sec = [intInf ; 180][2]
#
const lower_range = (x->sort(collect(x))[1]).([thickness_range, phase_range, buffer_range])
const upper_range = (x->sort(collect(x))[2]).([thickness_range, phase_range, buffer_range])
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
        ParticleSwarm(lower_range, upper_range, 3),          #1  Particle Swarm Optimization
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
println("")
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
global_objective_value = 1.0
global_count_stuck = 0
time_start_optim = 0.0
#
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
        CallbackInterval = 5,#1
        TargetFitness = 0.0, 
        MaxFuncEvals = max_steps,
        MaxTime = max_time_sec
    )
    #
    obj_result = best_fitness(optim_results)
    x_results  = best_candidate(optim_results)
    #
else
    #
    global time_start_optim = time()
    #
    function callback_optim(current_state)
        println(" * Efficiency = ", 1.0-current_state.value)
        flush(stdout)
        if current_state.value==global_objective_value
            global global_count_stuck+=1
        else
            global global_count_stuck=0
            global global_objective_value = current_state.value
        end
        if global_count_stuck>10 || time()-time_start_optim > max_time_sec
            return true
        else
            return false
        end 
    end
    #
    optim_results = Optim.optimize(objective_func,
        #lower_range, 
        #upper_range, 
        initial_guess_Optim,
        #Fminbox(chosen_solver), #This construction is bugged and gives the following error https://github.com/SciML/Optimization.jl/issues/179
        chosen_solver,
        Optim.Options(
            show_trace = true, 
            show_every = 1, 
            #If the callback returns true it stops the computation. If it returns false, it continues. 
            #It must return one of the two
            callback = callback_optim, 
            extended_trace = true,
            outer_iterations = max_steps,
            iterations = max_steps
            #,#f_calls_limit = 0
        )
    )
    #
    x_results = optim_results.minimizer
    obj_result = optim_results.minimum
    #
end
#
#Printing the final results
println("")
println("#"^25)
println("\n")        
println("** Optimized efficiency: " , 1.0-obj_result)
println("** Optimal parameters: ")
println("   1) disks_thickness     = "        , x_results[1])
println("   2) phase_shift         = "        , x_results[2])
println("   3) buffer_smooth       = "        , x_results[3])
#
#
@suppress begin global const r_lens = 6*lambda0 end
println("r_lens = ", r_lens)
println("** Optimal settings, but R=6: eta = ", SM_main(x_results[1], x_results[2], x_results[3]))
flush(stdout)



################## FUNCTIONS #############################################################################################################
#
function rand_range(range)
    sorted_range = sort(collect(range))
    (sorted_range[2]-sorted_range[1])*rand()+sorted_range[1]
end

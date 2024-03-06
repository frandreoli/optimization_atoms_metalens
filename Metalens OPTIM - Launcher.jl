#Designed for Julia higher than 1.5
using LinearAlgebra, Dates, HDF5, Random, Distributed
include("Metalens OPTIM - SM Functions.jl")
include("Metalens OPTIM - Core Evaluation.jl")
const intInf = typemax(Int)
#


################## OPTIONS ####################################################################################################################################
#
#System options
const z_fixed_option               =      ["YES" ; "NO"][1] #It must be 1
const phase_center_ring_option     =      ["YES" ; "NO"][1] #It must be 1
const fill_until_r_lens_option     =      [true ; false][1] #It must be 1
#
#Optimization options
const initial_guess_option         =      [true ; false][1] #If true, it feeds the algorithm with an initial guess
const monotonic_escape_option      =      [true ; false][1] #If true, it quits in case the optimization exhibits a non-monotonic behaviour
#
#Code generic options
const debug_option                 =      [nothing ; "ETA" ; "R_ATOMS"][2] 
#If "ETA", it computes the efficiency for an illustrative metalens and quits
#If "R_ATOMS", it saves the atomic positions for an illustrative metalens and quits

################## FIXED PARAMETERS ####################################################################################################################################
#
Random.seed!()
if nworkers()==1 addprocs() end
BLAS.set_num_threads(nworkers())
const TN = Float32
const k0 = 2.0*pi
const dipoles_polarization = [1.0 ; 0.0 ; 0.0]
const field_polarization   = [1.0 ; 0.0 ; 0.0]
#


################## GEOMETRICAL PARAMETERS ##############################################################################################################################
#
const w0                =    4.0
const focal_point       =    20.0
const r_lens            =    8.0
const gamma_prime       =    5.75
const laser_detuning    =    0.0
#
if debug_option == nothing
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
end
#


################## OPTIMIZATION PARAMETERS ################################################################################################
#
#Choice of the algorithm
solver_algorithm_index = 11
#
#From BlackBoxOptim (https://github.com/robertfeldt/BlackBoxOptim.jl):
#1  - Differential Evolution optimizer (metaheuristics) (suggested by developers)...(good)
#2  - Differential Evolution optimizer (metaheuristics) with limited radius.........(good)
#3  - Memetic algorithm.............................................................
#4  - Memetic algorithm with inheritance............................................(bad)
#5  - Stochastic Approximation algorithm............................................(worst)
#6  - Natural Evolution Strategies - separable......................................(good globally)
#7  - Natural Evolution Strategies - exponential....................................(very good globally)
#8  - Natural Evolution Strategies - exponential - distance-weighted................
#9  - Generating-set direct search..................................................(good in accuracy)
#10 - Generating-set direct search, with probabilistic descent.......................
#
#From Optim (https://julianlsolvers.github.io/Optim.jl/stable/):
#11 - Particle Swarm Optimization....................................................(best globally)
#12 - Simulated Annealing with intrinsic bounds......................................(not good: why?)
#13 - Simulated Annealing (bounds forced within the object function).................(not good)
#14 - Nelder-Mead (bounds forced within the object function).........................(not good)
#
#Generic parameters
const max_steps          = [intInf ; 2][1]
const max_time_sec       = [intInf ; 180][1]
const max_stuck_sequence = [intInf ; 50][1]
#
#Physical ranges and guess
const thickness_range    = (0.1,0.8)
const phase_range        = (-3.141592653589793,3.141592653589793)
const buffer_range       = (0.,0.5)
const initial_guess      = [0.6750965215553976, 1.1319090308946498, 0.46835704963894703] #[2.0/3, 1.0775,0.20048828125]


################## FUNCTIONS #############################################################################################################
#
function rand_range(range)
    sorted_range = sort(collect(range))
    (sorted_range[2]-sorted_range[1])*rand()+sorted_range[1]
end


################## OPTIMIZATION ###########################################################################################################
#
const lower_range = (x->sort(collect(x))[1]).([thickness_range, phase_range, buffer_range])
const upper_range = (x->sort(collect(x))[2]).([thickness_range, phase_range, buffer_range])
#
if debug_option == "R_ATOMS"
    debug_r_atoms(r_lens , focal_point, initial_guess[1], 0.8, initial_guess[2])
end
#
if debug_option == "ETA"
    parameters_debug     = [0.6750965215553976, 1.1319090308946498, 0.46835704963894703] #[0.704199570451851, 1.3429376990624875, 0.26280032985692664]
    debug_eta(10.0,w0 ,focal_point, gamma_prime , laser_detuning,parameters_debug)
end
#
#Definition of the solver
!(solver_algorithm_index in collect(1:14))  ? error("Undefined algorithm index") : nothing
solver_algorithm_index<=10 ? use_blackboxoptim_option = true : use_blackboxoptim_option = false
if use_blackboxoptim_option
    #
    using BlackBoxOptim
    chosen_solver = (
        :adaptive_de_rand_1_bin,                             
        :adaptive_de_rand_1_bin_radiuslimited,               
        :resampling_memetic_search,                           
        :resampling_inheritance_memetic_search,              
        :simultaneous_perturbation_stochastic_approximation, 
        :separable_nes,                                      
        :xnes,                                              
        :dxnes,                                                
        :generating_set_search,  
        :probabilistic_descent,                               
    )[solver_algorithm_index]
else
    using Optim
    chosen_solver = (
        ParticleSwarm(lower_range, upper_range, 70),          
        SAMIN(; rt=0.9),
        SimulatedAnnealing(),
        NelderMead()
    )[solver_algorithm_index-10]
end
#
#Printing data
println("")
println("The parameter ranges are: ")
println("thickness_range   = ", thickness_range)
println("phase_range       = ", phase_range)
println("buffer_range      = ", buffer_range)
if initial_guess_option
    println("The initial guesses are: ")
    println("- thickness guess = ", initial_guess[1])
    println("- phase guess     = "    , initial_guess[2])
    println("- buffer guess    = "   , initial_guess[3])
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
    return 1.0-SM_main(x[1], x[2], x[3], w0, focal_point, r_lens,gamma_prime,laser_detuning)
end
#
#Initializing the optimization
if initial_guess_option
    println("** STARTING the optimization WITH an initial guess **")
    blackboxoptim_input = (objective_func, initial_guess)
    initial_guess_Optim = initial_guess
else
    println("** STARTING the optimization WITHOUT an initial guess **")
    blackboxoptim_input = (objective_func,)
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
        CallbackFunction = x -> begin println(" * Efficiency = $(1-best_fitness(x)), $(best_candidate(x))\n"); flush(stdout) end,
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
        println(" * Efficiency = ", 1.0-current_state.value,"\n")
        flush(stdout)
        #
        if monotonic_escape_option
            if current_state.value>global_objective_value
                warn_string_monotonic = "The optimization has become non monotonic"
                println(("#"^25)*"\n****** WARNING: "*warn_string_monotonic*" *****\n"*("#"^25))
                @warn warn_string_monotonic
                return true
            end
        end
        #
        if current_state.value==global_objective_value
            global global_count_stuck+=1
        else
            global global_count_stuck=0
            global global_objective_value = current_state.value
        end
        if global_count_stuck>max_stuck_sequence || time()-time_start_optim > max_time_sec
            warn_string_stuck_or_time = "The optimization been stuck too many times with the same result OR time exceeded"
            println(("#"^25)*"\n****** WARNING: "*warn_string_stuck_or_time*" *****\n"*("#"^25))
            @warn warn_string_stuck_or_time
            return true
        else
            return false
        end 
    end
    #
    optim_options = Optim.Options(
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
    #
    if solver_algorithm_index<=12
        #
        optim_construct = (objective_func,lower_range,upper_range,initial_guess_Optim,chosen_solver,optim_options)
        #
    else
        #
        function objective_func_bounded(x)
            for index in 1:length(x)
                if x[index]>upper_range[index] || x[index]<lower_range[index]
                    return 10
                end
            end
            return objective_func(x)
        end
        #
        optim_construct = (objective_func_bounded,initial_guess_Optim,chosen_solver,optim_options)
    end
    #
    optim_results = Optim.optimize(optim_construct...)
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
println("** Optimized efficiency   = " , 1.0-obj_result)
println("** Optimal parameters: ")
println("   1) disks_thickness     = "        , x_results[1])
println("   2) phase_shift         = "        , x_results[2])
println("   3) buffer_smooth       = "        , x_results[3])
#
#
new_r_lens = 10.0#6.0
if r_lens<new_r_lens
    println("** Optimal settings, but r_lens = ",new_r_lens,": eta = ", SM_main(x_results[1], x_results[2], x_results[3], w0, focal_point, new_r_lens,gamma_prime,laser_detuning))
    flush(stdout)
end
 





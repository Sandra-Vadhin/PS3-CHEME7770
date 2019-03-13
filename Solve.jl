using GLPK
include("Flux.jl")
include("DataDictionary.jl")

# load the data dictionary -
data_dictionary = DataDictionary(0,100,1)

# solve the lp problem -
(objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag) = calculate_optimal_flux_distribution(data_dictionary)

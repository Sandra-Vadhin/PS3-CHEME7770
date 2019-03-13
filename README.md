# PS3-CHEME7770

Code requires GLPK to run. 
To check elemental balances of the cycle, run the following in the Julia REPL:

	julia> include("ElementalBalances.jl")
	
The FBA solution file Solve.jl calls Network.dat, DataDictionary.jl, and Flux.jl. To check the solution, enter the following in the Julia REPL:

	julia> include("Solve.jl")

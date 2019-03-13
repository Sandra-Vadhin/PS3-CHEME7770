# ----------------------------------------------------------------------------------- #
# Sandra Vadhin
# PS 3, CHEME 7770 Spring 2019 Cornell University
# March 12, 2019
# Adapted from JD Varner's code found at
#	https://github.com/varnerlab/CHEME7770-SimpleFBA-Problem
# ----------------------------------------------------------------------------------- #
# Copyright (c) 2017 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Function: DataDictionary
# Description: Holds simulation and model parameters as key => value pairs in a Julia Dict()
# Generated on: 2017-03-21T10:46:19.757
#
# Input arguments:
# time_start::Float64 => Simulation start time value (scalar)
# time_stop::Float64 => Simulation stop time value (scalar)
# time_step::Float64 => Simulation time step (scalar)
#
# Output arguments:
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model and simulation parameters as key => value pairs
# ----------------------------------------------------------------------------------- #
using DelimitedFiles

function DataDictionary(time_start,time_stop,time_step)

	# Load the stoichiometric network from disk -
	stoichiometric_matrix = readdlm("Network.dat");

	#calculate vmax values
	kcat_v1 = 203 * 3600;		#converting from s^-1 to hr^-1
	kcat_v2 = 34.5 * 3600;
	kcat_v3 = 249 * 3600;
	kcat_v4 = 88.1 * 3600;
	kcat_v5 = 13.7 * 3600;

	E = 0.01*10^-3;		#converting from umol/gdw to mmol/gdw

	vmax_1 = kcat_v1 * E;
	vmax_2 = kcat_v2 * E;
	vmax_3 = kcat_v3 * E;
	vmax_4 = kcat_v4 * E;
	vmax_5 = kcat_v5 * E;

	# Saturation function: [S]/(Km+[S]) when available
	sat_1 = (4.67E-3/(4.67E-3 + 3.92E-4))*(1.49E-2/(1.49E-2 + 1.54E-4));
	sat_2 = 1;
	sat_3 = 2.55E-4/(1.55E-3 + 2.55E-4)	;
	sat_4 = 4.49E-3/(1.60E-3 + 4.49E-3);
	sat_5 = 2.55E-4/(4.40E-03 + 2.55E-4);

	# Setup default flux bounds array -
	default_bounds_array = [

		0	vmax_1*sat_1;		# v1	ATP + L-citrulline + L-aspartate --> AMP + diphosphate + 2-(Nomega-L-arginino)succinate
		0	vmax_2*sat_2;		# v2 	2-(Nomega-L-arginino)succinate --> fumarate + L-arginine
		0	vmax_3*sat_3; # v3 	L-arginine + H2O --> L-ornithine + urea
		0	vmax_4*sat_4;		# v4 	carbamoyl phosphate + L-ornithine --> phosphate + L-citrulline
		-vmax_5	vmax_5*sat_5;		# v5 	2 L-arginine + 3 NADPH + 3 H+ + 4 O2 --> 2 L-citrulline + 2 nitric oxide + 3 NADP+ + 4 H2O
		0	10;		# b1 	input carbamoyl phosphate
		0	10;		# b2 	input aspartate
		0	10;		# b3  	output fumarate
		0	10;		# b4 	output urea
		0	10;		# b5 (input of ATP)
		0	10;		# b6 (output of AMP)
		0	10;		# b7 (output of diphosphate)
		-10	10;		# b8 (input of water)
		0	10;		# b9 (output of phosphate)
		-10	10;		# b10 (input of NADPH)
		-10	10;		# b11 (input of H+)
		-10	10;		# b12 (input of O2)
		-10	10;		# b13 nitric oxide out
		-10	10;		# b14 (output of NADP+)
		-10	10;		# b15 water out
		-10	10;		# b16 nitric oxide in
		];

	# Setup default species bounds array -
	species_bounds_array = [

		0.0	0.0	;	# 1 A_c
		0.0	0.0	;	# 2 A_x
		0.0	0.0	;	# 3 B_c
		0.0	0.0	;	# 4 B_x
		0.0	0.0	;	# 5 C_c
		0.0	0.0	;	# 6 C_x
		0.0	0.0	;	# 1 A_c
		0.0	0.0	;	# 2 A_x
		0.0	0.0	;	# 3 B_c
		0.0	0.0	;	# 4 B_x
		0.0	0.0	;	# 5 C_c
		0.0	0.0	;	# 6 C_x
		0.0	0.0	;	# 1 A_c
		0.0	0.0	;	# 2 A_x
		0.0	0.0	;	# 3 B_c
		0.0	0.0	;	# 4 B_x
		0.0	0.0	;	# 5 C_c
		0.0	0.0	;	# 6 C_x

	];

	# Setup the objective coefficient array -
	objective_coefficient_array = [

		0	;
		0	;
		0	;
		0	;
		0	;
		0	;
		0	;
		0	;
		1	;
		0	;
		0	;
		0	;
		0	;
		0	;
		0	;
		0	;
		0	;
		0	;
		0	;
	];

	# Min/Max flag - default is minimum -
	min_flag = false

	# List of reation strings - used to write flux report
	list_of_reaction_strings = [

	"v1: ATP + L-citrulline + L-aspartate --> AMP + diphosphate + 2-(Nomega-L-arginino)succinate"
	"v2: 2-(Nomega-L-arginino)succinate --> fumarate + L-arginine"
	"v3: L-arginine + H2O --> L-ornithine + urea"
	"v4: carbamoyl phosphate + L-ornithine --> phosphate + L-citrulline"
	"v5: 2 L-arginine + 3 NADPH + 3 H+ + 4 O2 --> 2 L-citrulline + 2 nitric oxide + 3 NADP+ + 4 H2O"
	"b1: input carbamoyl phosphate"
	"b2: input aspartate"
	"b3: output fumarate"
	"b4: output urea"
	"b5 (input of ATP)"
	"b6 (output of AMP)"
	"b7 (output of diphosphate)"
	"b8 (input of water)"
	"b9 (output of phosphate)"
	"b10 (input of NADPH)"
	"b11 (input of O2)"
	"b12 (output of H+)"
	"b13 (output of NO)"
	"b14 (output of NADP+)"
	"b15 water out"
	"b16 nitric oxide in"

	];

	# List of metabolite strings - used to write flux report
	list_of_metabolite_symbols = [

		"ATP"
		"citrulline"
		"aspartate"
		"AMP"
		"diphosphate"
		"argininosuccinate"
		"fumarate"
		"arginine"
		"water"
		"ornithine"
		"urea"
		"carbamoyl phosphate"
		"phosphate"
		"NADPH"
		"NADP+"
		"H+"
		"O2"
		"nitric oxide"
		"ADP"

	];

	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{AbstractString,Any}()
	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
	data_dictionary["objective_coefficient_array"] = objective_coefficient_array
	data_dictionary["default_flux_bounds_array"] = default_bounds_array;
	data_dictionary["species_bounds_array"] = species_bounds_array
	data_dictionary["list_of_reaction_strings"] = list_of_reaction_strings
	data_dictionary["list_of_metabolite_symbols"] = list_of_metabolite_symbols
	data_dictionary["min_flag"] = min_flag
	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #
	return data_dictionary
	return stoichiometric_matrix
	return objective_coefficient_array
	return default_bounds_array
	return species_bounds_array
	return min_flag
	return list_of_reaction_strings
	return list_of_metabolite_symbols
end

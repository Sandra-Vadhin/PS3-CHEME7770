#= ---------------------------------------------------------------------
Elemental Balance by Sandra Vadhin
for Problem Set 3, CHEME 7770 Spring 2019, Cornell University

Stoich rows 1-18, in order:
ATP, citrulline, aspartate, AMP, diphosphate,argininosuccinate, fumarate,
arginine, water, ornithine, urea, carbamoyl phosphate, phosphate,NADPH,
NADP+, H+, O2, nitric oxide

Stoich columns 1-21, in order:
v1, v2, v3, v4, v5, b1 (input of carbamoyl phosphate),
b2 (input of aspartate),b3 (output of fumarate), b4 (output of urea),
b5 (input of ATP), b6 (output of AMP), b7 (output of diphosphate),
b8 (input of water), b9 (output of phosphate), b10 (input of NADPH),
b11 (input of H+), b12 (output of O2), b13 (output of nitric oxide),
b14(output of NADP), b15 (output of H2O), b16 (input of nitric oxide)

Atoms rows 1-5: C,H,N,O,P

Atoms columns 1-18: compounds from the rows of Stoich
----------------------------------------------------------------------=#
using DelimitedFiles
Stoich = readdlm("Network.dat");

Atoms = [

10	6	4	10	0	10	4	6	0	5	1	1	0	21	21	0	0	0
16	13	7	14	4	18	4	14	2	12	4	4	3	30	29	1	0	0
5	3	1	5	0	4	0	4	0	2	2	1	0	7	7	0	0	1
13	3	4	7	7	6	4	2	1	2	1	5	4	17	17	0	2	1
3	0	0	1	2	0	0	0	0	0	0	1	1	3	3	0	0	0

]

Elements = Atoms*Stoich;
Rows = sum(Elements, dims = 2);

# v1-v5 are balanced. Overall, there is a deficit of one nitrogen.
# Rows gives a column vector for C,H,N,O,P

return Rows

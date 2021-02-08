## File Directory

SA.wrapper.r: R wrapper to calculate statistics (F_st, Pi_x, Pi_y) using ABC (approximate bayesian) 

Implicit Input parameter file: 'my.input.<n>'

Input file format:

the first is just n and is used by the ABC package, 
the recombination rate, 
the frequency of the SA allele on the X chromosome, 
the frequency of the SA allele on the Y, 
SAS location (in rho) from the SDR (in real world unknown)

Example Input files: 

A typical neutral simulation will therefore start with an input file (say input1) like:
1
100 
0
0
125

The only real difference with simulations with selection is that lines 3 and 4 will now have different numbers (let's say in input2 this time):
1
100
0.1
0.9
125

Note: Varying strengths of selection (so differences in frequency on X vs Y) as well as varying locations (final line) should be considered. For now let's not worry about recombination rate and keep that fixed on 100 rho/Mb.

The strength of SAS is determined by the third and fourth parameters of the simulation - the frequency of the SA site on the X and Y respectively. When these frequencies are very different (say 0.9 vs 0.1), then there is strong SAS selection. When they are both 0 there is none, and when they are equal to each other there is balancing selection (related to SAS, but distinct).

Explicit Input parameter: <n>

Explicit Output file: "outputn" with summary statistics

Usage: Rscript SA.wrapper.R <n>

Usage example: Rscript SA.wrapper.r 1


SaSimSnps: precompiled simulation binary code (called by the SA wrapper)

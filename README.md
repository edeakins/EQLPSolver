# EQLPSolver
EQLPSolver is the GitHub directory containing HiGHS 1-2-1 with Orbital Crossover implemented.

# Orbtial Crossover
Orbital Crossover is a method of obtaining a vertex solution to a Linear Programming instance from an interior point solution using symmetry information from the problem.

# Installation
Orbtial Crossover is implemented in HiGHS version 1.2.  With that, HiGHS requires cmake (version 3 mininum) for installation.  Once you have cmake installed the following steps will install HiGHS version 1.2 with orbital crossover capabilities on your machine.
* Clone the repository into a directory of your choice.
* cd into the EQLPSolver directory.
* Change the branch to the "stable" branch.  The default branch is "master".
* Change directory into "HiGHS-1-2-1".
* Make a directory named "build" and change directory into it.
* Run the command "cmake .." within the "build" directory.
* Once this command completes, run the command "make" within the "build" directory. 
* Once this command completes, the "highs" executable is located within "build/bin".

# Running the code
To run the code run the command "./highs filename" where filename is an ".lp" or ".mps" file.

# Running the code with solver options
The executable "./highs" can be given five different solver options via the command line argument "--solver". We list examples and explanations for the five solver options below.
* "./highs filename --solver=orbital_crossover_dual": Run HiGHS using orbital crossover with dual simplex to solve initial aggregate LP.
* "./highs filename --solver=orbital_crossover_ipm": Run HiGHS using orbital crossover with ipm and HiGHS crossover to solve initial aggregate LP.
* "./highs filename --solver=ipm_aggregate": Run HiGHS using ipm to solve initial aggregate LP, lift the solution to original LP, and use HiGHS crossover to achieve a vertex solution.
* "./highs filename --solver=dual": Run HiGHS using dual simplex solver to solve the LP.
* "./highs filename --solver=ipm": Run HiGHS using ipm and crossover to solve the LP.


# EQLPSolver
EQLPSolver is the GitHub directory containing HiGHS 1-2-1 with Orbital Crossover implemented.

# Orbtial Crossover
Orbital Crossover is a method of obtaining a vertex solution to a Linear Programming instance from an interior point solution using symmetry information from the problem.

# Installation
Orbtial Crossover is implemented in HiGHS version 1.2.  With that, HiGHS requires cmake (version 3 mininum) for installation.  EQLPSolver requires a unix operating system.  Once you have cmake installed the following steps will install HiGHS version 1.2 with orbital crossover capabilities on your machine.
* Clone the repository into a directory of your choice.
* cd into the EQLPSolver directory.
* Change the branch to the "stable" branch.  The default branch is "master".
* Change directory into "HiGHS-1-2-1".
* Make a directory named "build" and change directory into it.
  * Note: The "build" directory needs to be with a lower case "b" for the benchmarks scripts.
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

# Benchmark instances
Benchmark instances are stored in the directory "EQLPSolver/MIPLIB-2017" and "EQLPSolver/HS-COV-COD".  Initially, when you clone the repo, all instances are in a ".gz" file format.  You will need gunzip to decompress the files before passing them to the HiGHS executable.

# Running benchmarks
To benchmark and compare orbital crossover to other methods in HiGHS, we include shell scripts in the directory "EQLPSolver/TestScripts/" to be executed via a terminal.  In total, there are eleven scripts.  One of these scripts, "run_benchmarks.sh" will run HiGHS with all solver options for all instances in both the HS-COV-COD instance suite and the MIPLIB-2017 instance suite.  The other ten scripts will run an individual benchmark for a HiGHS solver option on one of the instance suites.  We list these below.
* "run_benchmarks_ocd_hs_cov_cod.sh": Run HiGHS with orbital crossover with dual simplex to solve initial aggregate LP on all instances in HS-COV-COD.
* "run_benchmarks_ocipm_hs_cov_cod.sh": Run HiGHS with orbital crossover with ipm and HiGHS crossover to solve initial aggregate LP on all instances in HS-COV-COD.
* "run_benchmarks_ipmagg_hs_cov_cod.sh": Run HiGHS with ipm and HiGHS crossover to solve initial aggregate LP, lift the solution to the original LP, and use HiGHS crossover to achieve a vertex solution on all instances in HS-COV_COD.
* "run_benchmarks_d_hs_cov_cod.sh": Run HiGHS with dual simplex on all instances in HS-COV-COD.
* "run_benchmarks_ipm_hs_cov_cod.sh": Run HiGHS with ipm and HiGHS crossover on all instances in HS-COV-COD.
* "run_benchmarks_ocd_miplib_2017.sh": Run HiGHS with orbital crossover with dual simplex to solve initial aggregate LP on all instances in MIPLIB-2017.
* "run_benchmarks_ocipm_miplib_2017.sh": Run HiGHS with orbital crossover with ipm and HiGHS crossover to solve initial aggregate LP on all instances in MIPLIB-2017.
* "run_benchmarks_ipmagg_miplib_2017.sh": Run HiGHS with ipm and HiGHS crossover to solve initial aggregate LP, lift the solution to the original LP, and use HiGHS crossover to achieve a vertex solution on all instances in MIPLIB-2017.
* "run_benchmarks_d_miplib_2017.sh": Run HiGHS with dual simplex on all instances in MIPLIB-2017.
* "run_benchmarks_ipm_miplib_2017.sh": Run HiGHS with ipm and HiGHS crossover on all instances in MIPLIB-2017.

All output files from any of these shell scripts will be placed in a directory named "Timings" inside whichever directory the shell script is called from.  For consistency, we recommend calling the shell scripts from the "EQLPSolver/TestSripts/" directory so that the output files will be contained in the "EQLPSolver/TestScripts/Timings/" directory.  Each output file is stored as ".csv" file and can be opened via Microsoft excel for easy viewing.  

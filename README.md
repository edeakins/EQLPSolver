# EQLPSolver
EQLPSolver is the GitHub directory containing HiGHS 1-2-1 with Orbital Crossover implemented.

# Orbtial Crossover
Orbital Crossover is a method of obtaining a vertex solution to a Linear Programming instance from an interior point solution using symmetry information from the problem.

# Installation
EQLPSolver requires a unix operating system.  Orbtial Crossover is implemented in HiGHS version 1.2.  HiGHS requires cmake (version 3 mininum) for installation.  Once you have cmake installed the following steps will install HiGHS version 1.2 with orbital crossover capabilities on your machine.
* Clone the repository into a directory of your choice.
* cd into the EQLPSolver directory.
* Change directory into "HiGHS-1-2-1".
* Make a directory named "build" and change directory into it.
  * Note: The "build" directory needs to be with a lower case "b" for the benchmarks scripts.
* Run the command "cmake .." within the "build" directory.
* Once this command completes, run the command "make" within the "build" directory. 
* Once this command completes, the "HiGHS" executable is located within "build/bin".

# Running the code
To run the code run the command "./highs filename" where filename is an ".lp" or ".mps" file.

# Running the code with solver options
The executable "./highs" can be given five different solver options via the command line argument "--solver". Two other command line options that might be frequently used are "--time_limit" that sets the time limit for the solver (in seconds) and "--report_time_file" that creates a "csv" file to record timings for specific components of the solver (specific to which "--solver" option is chosen). Below is an examle of how one might invoke the HiGHS executable with orbital crossover crossover as the solver, a timings file name "timings", and a time limit of one hour to solve a linear program (LP) stored in a file "lp.mps":
* ./highs "lp.mps" --solver=orbital_crossover_dual --report_time_file=timings --time_limit=3600  

# Benchmark instances
Benchmark instances are stored in the directory "EQLPSolver/MIPLIB-2017" and "EQLPSolver/HS-COV-COD".  Initially, when you clone the repo, all instances are in a ".gz" file format.  You will need gunzip to decompress the files before passing them to the HiGHS executable.

# Running benchmarks
To benchmark and compare orbital crossover to other methods in HiGHS, we include shell scripts in the directory "EQLPSolver/TestScripts/" to be executed via a terminal. There are many ".sh" files in this directory that exist so that individual solver methods can be ran for all instances in either the HS-COV-COD or MIPLIB-2017 directories. One ".sh" file named "run_benchmarks_all.sh" that will run the HiGHS executable for all instances in both the HS-COV-COD and MIPLIB-2017 directories for all solver methods. The other ".sh" files have a method number in the file name. This method number corresponds to a specific solver strategy. The mapping between method number, "--solver" option, and internal HiGHS solver strategy is below. Note that when using the ".sh" files to run the the test suite, you only need to run the ".sh" file itself, no command line options are necessary. However, we list the "--solver" command line options that one can use if calling the HiGHS executable below with the corresponding method number so that the user is aware of the possible inputs and what internal sovler strategy they trigger.

* Method 1 (command line -> --solver=simplex): Solve the given LP using HiGHS dual simplex implementation
* Method 2 (command line -> --solver=ipm): Solve the given LP using HiGHS interior point method (IPM) with HiGHS crossover
* Method 3 (command line -> --solver=orbital_crossover_dual): Solve the given LP using orbital crossover with iterative lifting; HiGHS dual simplex is used to solve the initial aggregate LP
* Method 4 (command line -> --solver=orbital_crossover_ipm): Solve the given LP using orbtial crossover with iterative lifting; HiGHS IPM with HiGHS crossover is used to solve the initial aggregate LP
* Method 5 (command line -> --solver=dual_alp_standard_crossover): Solve the given LP using HiGHS crossover with iterative lifting; HiGHS dual simplex is used to solve the initial aggregate LP
* Method 6 (command line -> --solver=dual_alp_primal_crossover): Solve the given LP using HiGHS primal crossover only with iterative lifting; HiGHS dual simplex is used to solve the initial aggregate LP
* Method 7 (command line -> --solver=ipmcross_alp_standard_crossover): Solve the given LP using HiGHS crossover with iterative lifting; HiGHS IPM and HiGHS crossover is used to solve the initial aggregate LP
* Method 8 (command line -> --solver=impcross_alp_primal_crossover): Solve the given LP using HiGHS primal crossover only with iterative lifting; HiGHS IPM and HiGHS crossover is used to solve the initial aggregate LP
* Method 9 (command line -> --solver=ipm_alp_standard_crossover): Solve the given LP using HiGHS crossover with iterative lifting; HiGHS IPM is used to solve the initial aggregate LP
* Method 10 (command line -> --solver=ipm_alp_primal_crossover): Solve the given LP using HiGHS primal crossover only with iterative lifting; HiGHS IPM is used to solve the initial aggregate LP
* Method 11 (command line -> --solver=orbital_crossover_dual_no_iter): Solve the given LP using orbital crossover without iterative lifting; HiGHS dual simplex is used to solve the initial aggregate LP
* Method 12 (command line -> --solver=orbital_crossover_ipm_no_iter): Solve the given LP using orbital crossover without iterative lifting; HiGHS IPM is used to solve the initial aggregate LP
* Method 13 (command line -> --solver=dual_alp_standard_crossover_iter): Solve the given LP using HiGHS crossover without iterative lifting; HiGHS dual simplex is used to solve the initial aggregate LP
* Method 14 (command line -> --solver=dual_alp_primal_crossover_iter): Solve the given LP using HiGHS primal crossover only without iterative lifting; HiGHS dual simplex is used to solve the initial aggregate LP
* Method 15 (command line -> --solver=ipmcross_alp_standard_crossover_iter): Solve the given LP using HiGHS crossover without iterative lifting; HiGHS IPM and HiGHS crossover is used to solve the initial aggregate LP
* Method 16 (command line -> --solver=ipmcross_alp_primal_crossover_iter): Solve the given LP using HiGHS primal crossover only without iterative lifting; HiGHS IPM and HiGHS crossover is used to solve the initial aggregate LP
* Method 17 (command line -> --solver=ipm_alp_standard_crossover_iter): Solve the given LP using HiGHS crossover without iterative lifting; HiGHS IPM is used to solve the initial aggregate LP
* Method 18 (command line -> --solver=ipm_alp_primal_crossover_iter): Solve the given LP using HiGHS primal crossover only without iterative lifting; HiGHS IPM is used to solve the initial aggregate LP



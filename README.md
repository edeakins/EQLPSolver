# Orbtial Crossover
Orbital Crossover is a method of obtaining a vertex solution to a Linear Programming instance from an interior point solution using symmetry information from the problem.

# Installation
Orbtial Crossover is implemented in HiGHS version 1.2.  With that, HiGHS requires cmake for installation.
* Clone the code into a directory of your choice.
* Change directory into "HiGHS-1-2-1".
* Make a directory named "build" and change directory into it.
* Run the command "cmake .." within the "build" directory.
* Once this command completes, run the command "make" within the "build" directory. 
* Once this command completes, the "highs" executable is located within "build/bin".

# Running the code
To run the code run the command "./highs filename --solver=orbital_crossover_dual" or "./highs filename --solver=orbital_crossover_dual".  The first method uses HiGHS serial dual simplex to solve ALP, and the second uses HiGHS interior point method with HiGHS crossover to solve ALP.

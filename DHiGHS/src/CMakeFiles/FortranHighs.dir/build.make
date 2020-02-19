# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.11

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.11.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.11.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS

# Include any dependencies generated for this target.
include src/CMakeFiles/FortranHighs.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/FortranHighs.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/FortranHighs.dir/flags.make

src/CMakeFiles/FortranHighs.dir/interfaces/highs_lp_solver.f90.o: src/CMakeFiles/FortranHighs.dir/flags.make
src/CMakeFiles/FortranHighs.dir/interfaces/highs_lp_solver.f90.o: src/interfaces/highs_lp_solver.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object src/CMakeFiles/FortranHighs.dir/interfaces/highs_lp_solver.f90.o"
	cd /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/src && /usr/local/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/src/interfaces/highs_lp_solver.f90 -o CMakeFiles/FortranHighs.dir/interfaces/highs_lp_solver.f90.o

src/CMakeFiles/FortranHighs.dir/interfaces/highs_lp_solver.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/FortranHighs.dir/interfaces/highs_lp_solver.f90.i"
	cd /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/src && /usr/local/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/src/interfaces/highs_lp_solver.f90 > CMakeFiles/FortranHighs.dir/interfaces/highs_lp_solver.f90.i

src/CMakeFiles/FortranHighs.dir/interfaces/highs_lp_solver.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/FortranHighs.dir/interfaces/highs_lp_solver.f90.s"
	cd /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/src && /usr/local/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/src/interfaces/highs_lp_solver.f90 -o CMakeFiles/FortranHighs.dir/interfaces/highs_lp_solver.f90.s

# Object files for target FortranHighs
FortranHighs_OBJECTS = \
"CMakeFiles/FortranHighs.dir/interfaces/highs_lp_solver.f90.o"

# External object files for target FortranHighs
FortranHighs_EXTERNAL_OBJECTS =

lib/libFortranHighs.dylib: src/CMakeFiles/FortranHighs.dir/interfaces/highs_lp_solver.f90.o
lib/libFortranHighs.dylib: src/CMakeFiles/FortranHighs.dir/build.make
lib/libFortranHighs.dylib: lib/libhighs.1.0.0.dylib
lib/libFortranHighs.dylib: src/CMakeFiles/FortranHighs.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran shared library ../lib/libFortranHighs.dylib"
	cd /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/FortranHighs.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/FortranHighs.dir/build: lib/libFortranHighs.dylib

.PHONY : src/CMakeFiles/FortranHighs.dir/build

src/CMakeFiles/FortranHighs.dir/clean:
	cd /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/src && $(CMAKE_COMMAND) -P CMakeFiles/FortranHighs.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/FortranHighs.dir/clean

src/CMakeFiles/FortranHighs.dir/depend:
	cd /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/src /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/src /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/src/CMakeFiles/FortranHighs.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/FortranHighs.dir/depend


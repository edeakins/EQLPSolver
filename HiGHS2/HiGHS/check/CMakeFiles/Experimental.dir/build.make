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

# Utility rule file for Experimental.

# Include the progress variables for this target.
include check/CMakeFiles/Experimental.dir/progress.make

check/CMakeFiles/Experimental:
	cd /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/check && /usr/local/Cellar/cmake/3.11.2/bin/ctest -D Experimental

Experimental: check/CMakeFiles/Experimental
Experimental: check/CMakeFiles/Experimental.dir/build.make

.PHONY : Experimental

# Rule to build all files generated by this target.
check/CMakeFiles/Experimental.dir/build: Experimental

.PHONY : check/CMakeFiles/Experimental.dir/build

check/CMakeFiles/Experimental.dir/clean:
	cd /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/check && $(CMAKE_COMMAND) -P CMakeFiles/Experimental.dir/cmake_clean.cmake
.PHONY : check/CMakeFiles/Experimental.dir/clean

check/CMakeFiles/Experimental.dir/depend:
	cd /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/check /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/check /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/check/CMakeFiles/Experimental.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : check/CMakeFiles/Experimental.dir/depend


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
CMAKE_SOURCE_DIR = /Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS/DevBuild

# Utility rule file for Continuous.

# Include the progress variables for this target.
include check/CMakeFiles/Continuous.dir/progress.make

check/CMakeFiles/Continuous:
	cd /Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/local/Cellar/cmake/3.11.2/bin/ctest -D Continuous

Continuous: check/CMakeFiles/Continuous
Continuous: check/CMakeFiles/Continuous.dir/build.make

.PHONY : Continuous

# Rule to build all files generated by this target.
check/CMakeFiles/Continuous.dir/build: Continuous

.PHONY : check/CMakeFiles/Continuous.dir/build

check/CMakeFiles/Continuous.dir/clean:
	cd /Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS/DevBuild/check && $(CMAKE_COMMAND) -P CMakeFiles/Continuous.dir/cmake_clean.cmake
.PHONY : check/CMakeFiles/Continuous.dir/clean

check/CMakeFiles/Continuous.dir/depend:
	cd /Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS/DevBuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS /Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS/check /Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS/DevBuild /Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS/DevBuild/check /Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS/DevBuild/check/CMakeFiles/Continuous.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : check/CMakeFiles/Continuous.dir/depend


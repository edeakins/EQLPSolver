# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/edeakins/LP/EQLPSolver/HiGHS-1-2-1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/edeakins/LP/EQLPSolver/HiGHS-1-2-1/debugBuild

# Utility rule file for Continuous.

# Include the progress variables for this target.
include check/CMakeFiles/Continuous.dir/progress.make

check/CMakeFiles/Continuous:
	cd /home/edeakins/LP/EQLPSolver/HiGHS-1-2-1/debugBuild/check && /usr/bin/ctest -D Continuous

Continuous: check/CMakeFiles/Continuous
Continuous: check/CMakeFiles/Continuous.dir/build.make

.PHONY : Continuous

# Rule to build all files generated by this target.
check/CMakeFiles/Continuous.dir/build: Continuous

.PHONY : check/CMakeFiles/Continuous.dir/build

check/CMakeFiles/Continuous.dir/clean:
	cd /home/edeakins/LP/EQLPSolver/HiGHS-1-2-1/debugBuild/check && $(CMAKE_COMMAND) -P CMakeFiles/Continuous.dir/cmake_clean.cmake
.PHONY : check/CMakeFiles/Continuous.dir/clean

check/CMakeFiles/Continuous.dir/depend:
	cd /home/edeakins/LP/EQLPSolver/HiGHS-1-2-1/debugBuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/edeakins/LP/EQLPSolver/HiGHS-1-2-1 /home/edeakins/LP/EQLPSolver/HiGHS-1-2-1/check /home/edeakins/LP/EQLPSolver/HiGHS-1-2-1/debugBuild /home/edeakins/LP/EQLPSolver/HiGHS-1-2-1/debugBuild/check /home/edeakins/LP/EQLPSolver/HiGHS-1-2-1/debugBuild/check/CMakeFiles/Continuous.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : check/CMakeFiles/Continuous.dir/depend


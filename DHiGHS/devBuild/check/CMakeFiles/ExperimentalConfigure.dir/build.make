# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.18.1/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.18.1/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS/devBuild

# Utility rule file for ExperimentalConfigure.

# Include the progress variables for this target.
include check/CMakeFiles/ExperimentalConfigure.dir/progress.make

check/CMakeFiles/ExperimentalConfigure:
	cd /Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS/devBuild/check && /usr/local/Cellar/cmake/3.18.1/bin/ctest -D ExperimentalConfigure

ExperimentalConfigure: check/CMakeFiles/ExperimentalConfigure
ExperimentalConfigure: check/CMakeFiles/ExperimentalConfigure.dir/build.make

.PHONY : ExperimentalConfigure

# Rule to build all files generated by this target.
check/CMakeFiles/ExperimentalConfigure.dir/build: ExperimentalConfigure

.PHONY : check/CMakeFiles/ExperimentalConfigure.dir/build

check/CMakeFiles/ExperimentalConfigure.dir/clean:
	cd /Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS/devBuild/check && $(CMAKE_COMMAND) -P CMakeFiles/ExperimentalConfigure.dir/cmake_clean.cmake
.PHONY : check/CMakeFiles/ExperimentalConfigure.dir/clean

check/CMakeFiles/ExperimentalConfigure.dir/depend:
	cd /Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS/devBuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS /Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS/check /Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS/devBuild /Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS/devBuild/check /Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS/devBuild/check/CMakeFiles/ExperimentalConfigure.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : check/CMakeFiles/ExperimentalConfigure.dir/depend


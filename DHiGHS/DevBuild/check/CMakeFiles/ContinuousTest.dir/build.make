# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_SOURCE_DIR = /home/edeakins/Thesis/EQLPSolver/DHiGHS

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild

# Utility rule file for ContinuousTest.

# Include the progress variables for this target.
include check/CMakeFiles/ContinuousTest.dir/progress.make

check/CMakeFiles/ContinuousTest:
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/ctest -D ContinuousTest

ContinuousTest: check/CMakeFiles/ContinuousTest
ContinuousTest: check/CMakeFiles/ContinuousTest.dir/build.make

.PHONY : ContinuousTest

# Rule to build all files generated by this target.
check/CMakeFiles/ContinuousTest.dir/build: ContinuousTest

.PHONY : check/CMakeFiles/ContinuousTest.dir/build

check/CMakeFiles/ContinuousTest.dir/clean:
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && $(CMAKE_COMMAND) -P CMakeFiles/ContinuousTest.dir/cmake_clean.cmake
.PHONY : check/CMakeFiles/ContinuousTest.dir/clean

check/CMakeFiles/ContinuousTest.dir/depend:
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/edeakins/Thesis/EQLPSolver/DHiGHS /home/edeakins/Thesis/EQLPSolver/DHiGHS/check /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check/CMakeFiles/ContinuousTest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : check/CMakeFiles/ContinuousTest.dir/depend


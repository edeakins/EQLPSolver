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
CMAKE_SOURCE_DIR = /home/edeakins/EQLPSolver/DHiGHS

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/edeakins/EQLPSolver/DHiGHS/DevBuild

# Utility rule file for NightlyTest.

# Include the progress variables for this target.
include check/CMakeFiles/NightlyTest.dir/progress.make

check/CMakeFiles/NightlyTest:
	cd /home/edeakins/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/ctest -D NightlyTest

NightlyTest: check/CMakeFiles/NightlyTest
NightlyTest: check/CMakeFiles/NightlyTest.dir/build.make

.PHONY : NightlyTest

# Rule to build all files generated by this target.
check/CMakeFiles/NightlyTest.dir/build: NightlyTest

.PHONY : check/CMakeFiles/NightlyTest.dir/build

check/CMakeFiles/NightlyTest.dir/clean:
	cd /home/edeakins/EQLPSolver/DHiGHS/DevBuild/check && $(CMAKE_COMMAND) -P CMakeFiles/NightlyTest.dir/cmake_clean.cmake
.PHONY : check/CMakeFiles/NightlyTest.dir/clean

check/CMakeFiles/NightlyTest.dir/depend:
	cd /home/edeakins/EQLPSolver/DHiGHS/DevBuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/edeakins/EQLPSolver/DHiGHS /home/edeakins/EQLPSolver/DHiGHS/check /home/edeakins/EQLPSolver/DHiGHS/DevBuild /home/edeakins/EQLPSolver/DHiGHS/DevBuild/check /home/edeakins/EQLPSolver/DHiGHS/DevBuild/check/CMakeFiles/NightlyTest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : check/CMakeFiles/NightlyTest.dir/depend


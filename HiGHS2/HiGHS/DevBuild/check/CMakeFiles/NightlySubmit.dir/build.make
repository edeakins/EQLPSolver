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
CMAKE_BINARY_DIR = /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/DevBuild

# Utility rule file for NightlySubmit.

# Include the progress variables for this target.
include check/CMakeFiles/NightlySubmit.dir/progress.make

check/CMakeFiles/NightlySubmit:
	cd /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/DevBuild/check && /usr/local/Cellar/cmake/3.11.2/bin/ctest -D NightlySubmit

NightlySubmit: check/CMakeFiles/NightlySubmit
NightlySubmit: check/CMakeFiles/NightlySubmit.dir/build.make

.PHONY : NightlySubmit

# Rule to build all files generated by this target.
check/CMakeFiles/NightlySubmit.dir/build: NightlySubmit

.PHONY : check/CMakeFiles/NightlySubmit.dir/build

check/CMakeFiles/NightlySubmit.dir/clean:
	cd /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/DevBuild/check && $(CMAKE_COMMAND) -P CMakeFiles/NightlySubmit.dir/cmake_clean.cmake
.PHONY : check/CMakeFiles/NightlySubmit.dir/clean

check/CMakeFiles/NightlySubmit.dir/depend:
	cd /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/DevBuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/check /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/DevBuild /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/DevBuild/check /Users/Chvatal/Documents/Thesis/EQLPSolver/HiGHS/DevBuild/check/CMakeFiles/NightlySubmit.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : check/CMakeFiles/NightlySubmit.dir/depend


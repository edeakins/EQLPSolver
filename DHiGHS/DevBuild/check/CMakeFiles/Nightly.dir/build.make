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
CMAKE_SOURCE_DIR = /home/edeak/EQLPSolver/DHiGHS

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/edeak/EQLPSolver/DHiGHS/DevBuild

# Utility rule file for Nightly.

# Include the progress variables for this target.
include check/CMakeFiles/Nightly.dir/progress.make

check/CMakeFiles/Nightly:
	cd /home/edeak/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/ctest -D Nightly

Nightly: check/CMakeFiles/Nightly
Nightly: check/CMakeFiles/Nightly.dir/build.make

.PHONY : Nightly

# Rule to build all files generated by this target.
check/CMakeFiles/Nightly.dir/build: Nightly

.PHONY : check/CMakeFiles/Nightly.dir/build

check/CMakeFiles/Nightly.dir/clean:
	cd /home/edeak/EQLPSolver/DHiGHS/DevBuild/check && $(CMAKE_COMMAND) -P CMakeFiles/Nightly.dir/cmake_clean.cmake
.PHONY : check/CMakeFiles/Nightly.dir/clean

check/CMakeFiles/Nightly.dir/depend:
	cd /home/edeak/EQLPSolver/DHiGHS/DevBuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/edeak/EQLPSolver/DHiGHS /home/edeak/EQLPSolver/DHiGHS/check /home/edeak/EQLPSolver/DHiGHS/DevBuild /home/edeak/EQLPSolver/DHiGHS/DevBuild/check /home/edeak/EQLPSolver/DHiGHS/DevBuild/check/CMakeFiles/Nightly.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : check/CMakeFiles/Nightly.dir/depend


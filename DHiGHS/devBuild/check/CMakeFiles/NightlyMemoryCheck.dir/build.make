# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

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
CMAKE_COMMAND = /Applications/CMake.app/Contents/bin/cmake

# The command to remove a file.
RM = /Applications/CMake.app/Contents/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/ethanjedidahdeakins/Work/LP/EQLPSolver/DHiGHS

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/ethanjedidahdeakins/Work/LP/EQLPSolver/DHiGHS/devBuild

# Utility rule file for NightlyMemoryCheck.

# Include any custom commands dependencies for this target.
include check/CMakeFiles/NightlyMemoryCheck.dir/compiler_depend.make

# Include the progress variables for this target.
include check/CMakeFiles/NightlyMemoryCheck.dir/progress.make

check/CMakeFiles/NightlyMemoryCheck:
	cd /Users/ethanjedidahdeakins/Work/LP/EQLPSolver/DHiGHS/devBuild/check && /Applications/CMake.app/Contents/bin/ctest -D NightlyMemoryCheck

NightlyMemoryCheck: check/CMakeFiles/NightlyMemoryCheck
NightlyMemoryCheck: check/CMakeFiles/NightlyMemoryCheck.dir/build.make
.PHONY : NightlyMemoryCheck

# Rule to build all files generated by this target.
check/CMakeFiles/NightlyMemoryCheck.dir/build: NightlyMemoryCheck
.PHONY : check/CMakeFiles/NightlyMemoryCheck.dir/build

check/CMakeFiles/NightlyMemoryCheck.dir/clean:
	cd /Users/ethanjedidahdeakins/Work/LP/EQLPSolver/DHiGHS/devBuild/check && $(CMAKE_COMMAND) -P CMakeFiles/NightlyMemoryCheck.dir/cmake_clean.cmake
.PHONY : check/CMakeFiles/NightlyMemoryCheck.dir/clean

check/CMakeFiles/NightlyMemoryCheck.dir/depend:
	cd /Users/ethanjedidahdeakins/Work/LP/EQLPSolver/DHiGHS/devBuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ethanjedidahdeakins/Work/LP/EQLPSolver/DHiGHS /Users/ethanjedidahdeakins/Work/LP/EQLPSolver/DHiGHS/check /Users/ethanjedidahdeakins/Work/LP/EQLPSolver/DHiGHS/devBuild /Users/ethanjedidahdeakins/Work/LP/EQLPSolver/DHiGHS/devBuild/check /Users/ethanjedidahdeakins/Work/LP/EQLPSolver/DHiGHS/devBuild/check/CMakeFiles/NightlyMemoryCheck.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : check/CMakeFiles/NightlyMemoryCheck.dir/depend


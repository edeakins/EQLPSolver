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
CMAKE_SOURCE_DIR = /home/edeakins/LP/EQLPSolver/DHiGHS

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/edeakins/LP/EQLPSolver/DHiGHS/debugBuild

# Include any dependencies generated for this target.
include check/CMakeFiles/capi_unit_tests.dir/depend.make

# Include the progress variables for this target.
include check/CMakeFiles/capi_unit_tests.dir/progress.make

# Include the compile flags for this target's objects.
include check/CMakeFiles/capi_unit_tests.dir/flags.make

check/CMakeFiles/capi_unit_tests.dir/TestCAPI.c.o: check/CMakeFiles/capi_unit_tests.dir/flags.make
check/CMakeFiles/capi_unit_tests.dir/TestCAPI.c.o: ../check/TestCAPI.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/edeakins/LP/EQLPSolver/DHiGHS/debugBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object check/CMakeFiles/capi_unit_tests.dir/TestCAPI.c.o"
	cd /home/edeakins/LP/EQLPSolver/DHiGHS/debugBuild/check && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/capi_unit_tests.dir/TestCAPI.c.o   -c /home/edeakins/LP/EQLPSolver/DHiGHS/check/TestCAPI.c

check/CMakeFiles/capi_unit_tests.dir/TestCAPI.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/capi_unit_tests.dir/TestCAPI.c.i"
	cd /home/edeakins/LP/EQLPSolver/DHiGHS/debugBuild/check && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/edeakins/LP/EQLPSolver/DHiGHS/check/TestCAPI.c > CMakeFiles/capi_unit_tests.dir/TestCAPI.c.i

check/CMakeFiles/capi_unit_tests.dir/TestCAPI.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/capi_unit_tests.dir/TestCAPI.c.s"
	cd /home/edeakins/LP/EQLPSolver/DHiGHS/debugBuild/check && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/edeakins/LP/EQLPSolver/DHiGHS/check/TestCAPI.c -o CMakeFiles/capi_unit_tests.dir/TestCAPI.c.s

# Object files for target capi_unit_tests
capi_unit_tests_OBJECTS = \
"CMakeFiles/capi_unit_tests.dir/TestCAPI.c.o"

# External object files for target capi_unit_tests
capi_unit_tests_EXTERNAL_OBJECTS =

bin/capi_unit_tests: check/CMakeFiles/capi_unit_tests.dir/TestCAPI.c.o
bin/capi_unit_tests: check/CMakeFiles/capi_unit_tests.dir/build.make
bin/capi_unit_tests: lib/libhighs.so.1.0.0
bin/capi_unit_tests: lib/libipx.so
bin/capi_unit_tests: lib/libbasiclu.so
bin/capi_unit_tests: check/CMakeFiles/capi_unit_tests.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/edeakins/LP/EQLPSolver/DHiGHS/debugBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable ../bin/capi_unit_tests"
	cd /home/edeakins/LP/EQLPSolver/DHiGHS/debugBuild/check && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/capi_unit_tests.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
check/CMakeFiles/capi_unit_tests.dir/build: bin/capi_unit_tests

.PHONY : check/CMakeFiles/capi_unit_tests.dir/build

check/CMakeFiles/capi_unit_tests.dir/clean:
	cd /home/edeakins/LP/EQLPSolver/DHiGHS/debugBuild/check && $(CMAKE_COMMAND) -P CMakeFiles/capi_unit_tests.dir/cmake_clean.cmake
.PHONY : check/CMakeFiles/capi_unit_tests.dir/clean

check/CMakeFiles/capi_unit_tests.dir/depend:
	cd /home/edeakins/LP/EQLPSolver/DHiGHS/debugBuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/edeakins/LP/EQLPSolver/DHiGHS /home/edeakins/LP/EQLPSolver/DHiGHS/check /home/edeakins/LP/EQLPSolver/DHiGHS/debugBuild /home/edeakins/LP/EQLPSolver/DHiGHS/debugBuild/check /home/edeakins/LP/EQLPSolver/DHiGHS/debugBuild/check/CMakeFiles/capi_unit_tests.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : check/CMakeFiles/capi_unit_tests.dir/depend

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

# Include any dependencies generated for this target.
include check/CMakeFiles/unit_tests.dir/depend.make

# Include the progress variables for this target.
include check/CMakeFiles/unit_tests.dir/progress.make

# Include the compile flags for this target's objects.
include check/CMakeFiles/unit_tests.dir/flags.make

check/CMakeFiles/unit_tests.dir/TestMain.cpp.o: check/CMakeFiles/unit_tests.dir/flags.make
check/CMakeFiles/unit_tests.dir/TestMain.cpp.o: ../check/TestMain.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object check/CMakeFiles/unit_tests.dir/TestMain.cpp.o"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/unit_tests.dir/TestMain.cpp.o -c /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestMain.cpp

check/CMakeFiles/unit_tests.dir/TestMain.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/unit_tests.dir/TestMain.cpp.i"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestMain.cpp > CMakeFiles/unit_tests.dir/TestMain.cpp.i

check/CMakeFiles/unit_tests.dir/TestMain.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/unit_tests.dir/TestMain.cpp.s"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestMain.cpp -o CMakeFiles/unit_tests.dir/TestMain.cpp.s

check/CMakeFiles/unit_tests.dir/TestMain.cpp.o.requires:

.PHONY : check/CMakeFiles/unit_tests.dir/TestMain.cpp.o.requires

check/CMakeFiles/unit_tests.dir/TestMain.cpp.o.provides: check/CMakeFiles/unit_tests.dir/TestMain.cpp.o.requires
	$(MAKE) -f check/CMakeFiles/unit_tests.dir/build.make check/CMakeFiles/unit_tests.dir/TestMain.cpp.o.provides.build
.PHONY : check/CMakeFiles/unit_tests.dir/TestMain.cpp.o.provides

check/CMakeFiles/unit_tests.dir/TestMain.cpp.o.provides.build: check/CMakeFiles/unit_tests.dir/TestMain.cpp.o


check/CMakeFiles/unit_tests.dir/TestOptions.cpp.o: check/CMakeFiles/unit_tests.dir/flags.make
check/CMakeFiles/unit_tests.dir/TestOptions.cpp.o: ../check/TestOptions.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object check/CMakeFiles/unit_tests.dir/TestOptions.cpp.o"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/unit_tests.dir/TestOptions.cpp.o -c /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestOptions.cpp

check/CMakeFiles/unit_tests.dir/TestOptions.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/unit_tests.dir/TestOptions.cpp.i"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestOptions.cpp > CMakeFiles/unit_tests.dir/TestOptions.cpp.i

check/CMakeFiles/unit_tests.dir/TestOptions.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/unit_tests.dir/TestOptions.cpp.s"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestOptions.cpp -o CMakeFiles/unit_tests.dir/TestOptions.cpp.s

check/CMakeFiles/unit_tests.dir/TestOptions.cpp.o.requires:

.PHONY : check/CMakeFiles/unit_tests.dir/TestOptions.cpp.o.requires

check/CMakeFiles/unit_tests.dir/TestOptions.cpp.o.provides: check/CMakeFiles/unit_tests.dir/TestOptions.cpp.o.requires
	$(MAKE) -f check/CMakeFiles/unit_tests.dir/build.make check/CMakeFiles/unit_tests.dir/TestOptions.cpp.o.provides.build
.PHONY : check/CMakeFiles/unit_tests.dir/TestOptions.cpp.o.provides

check/CMakeFiles/unit_tests.dir/TestOptions.cpp.o.provides.build: check/CMakeFiles/unit_tests.dir/TestOptions.cpp.o


check/CMakeFiles/unit_tests.dir/TestPresolve.cpp.o: check/CMakeFiles/unit_tests.dir/flags.make
check/CMakeFiles/unit_tests.dir/TestPresolve.cpp.o: ../check/TestPresolve.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object check/CMakeFiles/unit_tests.dir/TestPresolve.cpp.o"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/unit_tests.dir/TestPresolve.cpp.o -c /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestPresolve.cpp

check/CMakeFiles/unit_tests.dir/TestPresolve.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/unit_tests.dir/TestPresolve.cpp.i"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestPresolve.cpp > CMakeFiles/unit_tests.dir/TestPresolve.cpp.i

check/CMakeFiles/unit_tests.dir/TestPresolve.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/unit_tests.dir/TestPresolve.cpp.s"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestPresolve.cpp -o CMakeFiles/unit_tests.dir/TestPresolve.cpp.s

check/CMakeFiles/unit_tests.dir/TestPresolve.cpp.o.requires:

.PHONY : check/CMakeFiles/unit_tests.dir/TestPresolve.cpp.o.requires

check/CMakeFiles/unit_tests.dir/TestPresolve.cpp.o.provides: check/CMakeFiles/unit_tests.dir/TestPresolve.cpp.o.requires
	$(MAKE) -f check/CMakeFiles/unit_tests.dir/build.make check/CMakeFiles/unit_tests.dir/TestPresolve.cpp.o.provides.build
.PHONY : check/CMakeFiles/unit_tests.dir/TestPresolve.cpp.o.provides

check/CMakeFiles/unit_tests.dir/TestPresolve.cpp.o.provides.build: check/CMakeFiles/unit_tests.dir/TestPresolve.cpp.o


check/CMakeFiles/unit_tests.dir/TestIO.cpp.o: check/CMakeFiles/unit_tests.dir/flags.make
check/CMakeFiles/unit_tests.dir/TestIO.cpp.o: ../check/TestIO.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object check/CMakeFiles/unit_tests.dir/TestIO.cpp.o"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/unit_tests.dir/TestIO.cpp.o -c /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestIO.cpp

check/CMakeFiles/unit_tests.dir/TestIO.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/unit_tests.dir/TestIO.cpp.i"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestIO.cpp > CMakeFiles/unit_tests.dir/TestIO.cpp.i

check/CMakeFiles/unit_tests.dir/TestIO.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/unit_tests.dir/TestIO.cpp.s"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestIO.cpp -o CMakeFiles/unit_tests.dir/TestIO.cpp.s

check/CMakeFiles/unit_tests.dir/TestIO.cpp.o.requires:

.PHONY : check/CMakeFiles/unit_tests.dir/TestIO.cpp.o.requires

check/CMakeFiles/unit_tests.dir/TestIO.cpp.o.provides: check/CMakeFiles/unit_tests.dir/TestIO.cpp.o.requires
	$(MAKE) -f check/CMakeFiles/unit_tests.dir/build.make check/CMakeFiles/unit_tests.dir/TestIO.cpp.o.provides.build
.PHONY : check/CMakeFiles/unit_tests.dir/TestIO.cpp.o.provides

check/CMakeFiles/unit_tests.dir/TestIO.cpp.o.provides.build: check/CMakeFiles/unit_tests.dir/TestIO.cpp.o


check/CMakeFiles/unit_tests.dir/TestSort.cpp.o: check/CMakeFiles/unit_tests.dir/flags.make
check/CMakeFiles/unit_tests.dir/TestSort.cpp.o: ../check/TestSort.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object check/CMakeFiles/unit_tests.dir/TestSort.cpp.o"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/unit_tests.dir/TestSort.cpp.o -c /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestSort.cpp

check/CMakeFiles/unit_tests.dir/TestSort.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/unit_tests.dir/TestSort.cpp.i"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestSort.cpp > CMakeFiles/unit_tests.dir/TestSort.cpp.i

check/CMakeFiles/unit_tests.dir/TestSort.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/unit_tests.dir/TestSort.cpp.s"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestSort.cpp -o CMakeFiles/unit_tests.dir/TestSort.cpp.s

check/CMakeFiles/unit_tests.dir/TestSort.cpp.o.requires:

.PHONY : check/CMakeFiles/unit_tests.dir/TestSort.cpp.o.requires

check/CMakeFiles/unit_tests.dir/TestSort.cpp.o.provides: check/CMakeFiles/unit_tests.dir/TestSort.cpp.o.requires
	$(MAKE) -f check/CMakeFiles/unit_tests.dir/build.make check/CMakeFiles/unit_tests.dir/TestSort.cpp.o.provides.build
.PHONY : check/CMakeFiles/unit_tests.dir/TestSort.cpp.o.provides

check/CMakeFiles/unit_tests.dir/TestSort.cpp.o.provides.build: check/CMakeFiles/unit_tests.dir/TestSort.cpp.o


check/CMakeFiles/unit_tests.dir/TestSetup.cpp.o: check/CMakeFiles/unit_tests.dir/flags.make
check/CMakeFiles/unit_tests.dir/TestSetup.cpp.o: ../check/TestSetup.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object check/CMakeFiles/unit_tests.dir/TestSetup.cpp.o"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/unit_tests.dir/TestSetup.cpp.o -c /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestSetup.cpp

check/CMakeFiles/unit_tests.dir/TestSetup.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/unit_tests.dir/TestSetup.cpp.i"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestSetup.cpp > CMakeFiles/unit_tests.dir/TestSetup.cpp.i

check/CMakeFiles/unit_tests.dir/TestSetup.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/unit_tests.dir/TestSetup.cpp.s"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestSetup.cpp -o CMakeFiles/unit_tests.dir/TestSetup.cpp.s

check/CMakeFiles/unit_tests.dir/TestSetup.cpp.o.requires:

.PHONY : check/CMakeFiles/unit_tests.dir/TestSetup.cpp.o.requires

check/CMakeFiles/unit_tests.dir/TestSetup.cpp.o.provides: check/CMakeFiles/unit_tests.dir/TestSetup.cpp.o.requires
	$(MAKE) -f check/CMakeFiles/unit_tests.dir/build.make check/CMakeFiles/unit_tests.dir/TestSetup.cpp.o.provides.build
.PHONY : check/CMakeFiles/unit_tests.dir/TestSetup.cpp.o.provides

check/CMakeFiles/unit_tests.dir/TestSetup.cpp.o.provides.build: check/CMakeFiles/unit_tests.dir/TestSetup.cpp.o


check/CMakeFiles/unit_tests.dir/TestFilereader.cpp.o: check/CMakeFiles/unit_tests.dir/flags.make
check/CMakeFiles/unit_tests.dir/TestFilereader.cpp.o: ../check/TestFilereader.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object check/CMakeFiles/unit_tests.dir/TestFilereader.cpp.o"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/unit_tests.dir/TestFilereader.cpp.o -c /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestFilereader.cpp

check/CMakeFiles/unit_tests.dir/TestFilereader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/unit_tests.dir/TestFilereader.cpp.i"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestFilereader.cpp > CMakeFiles/unit_tests.dir/TestFilereader.cpp.i

check/CMakeFiles/unit_tests.dir/TestFilereader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/unit_tests.dir/TestFilereader.cpp.s"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestFilereader.cpp -o CMakeFiles/unit_tests.dir/TestFilereader.cpp.s

check/CMakeFiles/unit_tests.dir/TestFilereader.cpp.o.requires:

.PHONY : check/CMakeFiles/unit_tests.dir/TestFilereader.cpp.o.requires

check/CMakeFiles/unit_tests.dir/TestFilereader.cpp.o.provides: check/CMakeFiles/unit_tests.dir/TestFilereader.cpp.o.requires
	$(MAKE) -f check/CMakeFiles/unit_tests.dir/build.make check/CMakeFiles/unit_tests.dir/TestFilereader.cpp.o.provides.build
.PHONY : check/CMakeFiles/unit_tests.dir/TestFilereader.cpp.o.provides

check/CMakeFiles/unit_tests.dir/TestFilereader.cpp.o.provides.build: check/CMakeFiles/unit_tests.dir/TestFilereader.cpp.o


check/CMakeFiles/unit_tests.dir/TestInfo.cpp.o: check/CMakeFiles/unit_tests.dir/flags.make
check/CMakeFiles/unit_tests.dir/TestInfo.cpp.o: ../check/TestInfo.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object check/CMakeFiles/unit_tests.dir/TestInfo.cpp.o"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/unit_tests.dir/TestInfo.cpp.o -c /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestInfo.cpp

check/CMakeFiles/unit_tests.dir/TestInfo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/unit_tests.dir/TestInfo.cpp.i"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestInfo.cpp > CMakeFiles/unit_tests.dir/TestInfo.cpp.i

check/CMakeFiles/unit_tests.dir/TestInfo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/unit_tests.dir/TestInfo.cpp.s"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestInfo.cpp -o CMakeFiles/unit_tests.dir/TestInfo.cpp.s

check/CMakeFiles/unit_tests.dir/TestInfo.cpp.o.requires:

.PHONY : check/CMakeFiles/unit_tests.dir/TestInfo.cpp.o.requires

check/CMakeFiles/unit_tests.dir/TestInfo.cpp.o.provides: check/CMakeFiles/unit_tests.dir/TestInfo.cpp.o.requires
	$(MAKE) -f check/CMakeFiles/unit_tests.dir/build.make check/CMakeFiles/unit_tests.dir/TestInfo.cpp.o.provides.build
.PHONY : check/CMakeFiles/unit_tests.dir/TestInfo.cpp.o.provides

check/CMakeFiles/unit_tests.dir/TestInfo.cpp.o.provides.build: check/CMakeFiles/unit_tests.dir/TestInfo.cpp.o


check/CMakeFiles/unit_tests.dir/TestBasisSolves.cpp.o: check/CMakeFiles/unit_tests.dir/flags.make
check/CMakeFiles/unit_tests.dir/TestBasisSolves.cpp.o: ../check/TestBasisSolves.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object check/CMakeFiles/unit_tests.dir/TestBasisSolves.cpp.o"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/unit_tests.dir/TestBasisSolves.cpp.o -c /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestBasisSolves.cpp

check/CMakeFiles/unit_tests.dir/TestBasisSolves.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/unit_tests.dir/TestBasisSolves.cpp.i"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestBasisSolves.cpp > CMakeFiles/unit_tests.dir/TestBasisSolves.cpp.i

check/CMakeFiles/unit_tests.dir/TestBasisSolves.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/unit_tests.dir/TestBasisSolves.cpp.s"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestBasisSolves.cpp -o CMakeFiles/unit_tests.dir/TestBasisSolves.cpp.s

check/CMakeFiles/unit_tests.dir/TestBasisSolves.cpp.o.requires:

.PHONY : check/CMakeFiles/unit_tests.dir/TestBasisSolves.cpp.o.requires

check/CMakeFiles/unit_tests.dir/TestBasisSolves.cpp.o.provides: check/CMakeFiles/unit_tests.dir/TestBasisSolves.cpp.o.requires
	$(MAKE) -f check/CMakeFiles/unit_tests.dir/build.make check/CMakeFiles/unit_tests.dir/TestBasisSolves.cpp.o.provides.build
.PHONY : check/CMakeFiles/unit_tests.dir/TestBasisSolves.cpp.o.provides

check/CMakeFiles/unit_tests.dir/TestBasisSolves.cpp.o.provides.build: check/CMakeFiles/unit_tests.dir/TestBasisSolves.cpp.o


check/CMakeFiles/unit_tests.dir/TestLpValidation.cpp.o: check/CMakeFiles/unit_tests.dir/flags.make
check/CMakeFiles/unit_tests.dir/TestLpValidation.cpp.o: ../check/TestLpValidation.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object check/CMakeFiles/unit_tests.dir/TestLpValidation.cpp.o"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/unit_tests.dir/TestLpValidation.cpp.o -c /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestLpValidation.cpp

check/CMakeFiles/unit_tests.dir/TestLpValidation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/unit_tests.dir/TestLpValidation.cpp.i"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestLpValidation.cpp > CMakeFiles/unit_tests.dir/TestLpValidation.cpp.i

check/CMakeFiles/unit_tests.dir/TestLpValidation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/unit_tests.dir/TestLpValidation.cpp.s"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestLpValidation.cpp -o CMakeFiles/unit_tests.dir/TestLpValidation.cpp.s

check/CMakeFiles/unit_tests.dir/TestLpValidation.cpp.o.requires:

.PHONY : check/CMakeFiles/unit_tests.dir/TestLpValidation.cpp.o.requires

check/CMakeFiles/unit_tests.dir/TestLpValidation.cpp.o.provides: check/CMakeFiles/unit_tests.dir/TestLpValidation.cpp.o.requires
	$(MAKE) -f check/CMakeFiles/unit_tests.dir/build.make check/CMakeFiles/unit_tests.dir/TestLpValidation.cpp.o.provides.build
.PHONY : check/CMakeFiles/unit_tests.dir/TestLpValidation.cpp.o.provides

check/CMakeFiles/unit_tests.dir/TestLpValidation.cpp.o.provides.build: check/CMakeFiles/unit_tests.dir/TestLpValidation.cpp.o


check/CMakeFiles/unit_tests.dir/TestLpModification.cpp.o: check/CMakeFiles/unit_tests.dir/flags.make
check/CMakeFiles/unit_tests.dir/TestLpModification.cpp.o: ../check/TestLpModification.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object check/CMakeFiles/unit_tests.dir/TestLpModification.cpp.o"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/unit_tests.dir/TestLpModification.cpp.o -c /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestLpModification.cpp

check/CMakeFiles/unit_tests.dir/TestLpModification.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/unit_tests.dir/TestLpModification.cpp.i"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestLpModification.cpp > CMakeFiles/unit_tests.dir/TestLpModification.cpp.i

check/CMakeFiles/unit_tests.dir/TestLpModification.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/unit_tests.dir/TestLpModification.cpp.s"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestLpModification.cpp -o CMakeFiles/unit_tests.dir/TestLpModification.cpp.s

check/CMakeFiles/unit_tests.dir/TestLpModification.cpp.o.requires:

.PHONY : check/CMakeFiles/unit_tests.dir/TestLpModification.cpp.o.requires

check/CMakeFiles/unit_tests.dir/TestLpModification.cpp.o.provides: check/CMakeFiles/unit_tests.dir/TestLpModification.cpp.o.requires
	$(MAKE) -f check/CMakeFiles/unit_tests.dir/build.make check/CMakeFiles/unit_tests.dir/TestLpModification.cpp.o.provides.build
.PHONY : check/CMakeFiles/unit_tests.dir/TestLpModification.cpp.o.provides

check/CMakeFiles/unit_tests.dir/TestLpModification.cpp.o.provides.build: check/CMakeFiles/unit_tests.dir/TestLpModification.cpp.o


check/CMakeFiles/unit_tests.dir/Avgas.cpp.o: check/CMakeFiles/unit_tests.dir/flags.make
check/CMakeFiles/unit_tests.dir/Avgas.cpp.o: ../check/Avgas.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object check/CMakeFiles/unit_tests.dir/Avgas.cpp.o"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/unit_tests.dir/Avgas.cpp.o -c /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/Avgas.cpp

check/CMakeFiles/unit_tests.dir/Avgas.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/unit_tests.dir/Avgas.cpp.i"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/Avgas.cpp > CMakeFiles/unit_tests.dir/Avgas.cpp.i

check/CMakeFiles/unit_tests.dir/Avgas.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/unit_tests.dir/Avgas.cpp.s"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/Avgas.cpp -o CMakeFiles/unit_tests.dir/Avgas.cpp.s

check/CMakeFiles/unit_tests.dir/Avgas.cpp.o.requires:

.PHONY : check/CMakeFiles/unit_tests.dir/Avgas.cpp.o.requires

check/CMakeFiles/unit_tests.dir/Avgas.cpp.o.provides: check/CMakeFiles/unit_tests.dir/Avgas.cpp.o.requires
	$(MAKE) -f check/CMakeFiles/unit_tests.dir/build.make check/CMakeFiles/unit_tests.dir/Avgas.cpp.o.provides.build
.PHONY : check/CMakeFiles/unit_tests.dir/Avgas.cpp.o.provides

check/CMakeFiles/unit_tests.dir/Avgas.cpp.o.provides.build: check/CMakeFiles/unit_tests.dir/Avgas.cpp.o


check/CMakeFiles/unit_tests.dir/TestIpx.cpp.o: check/CMakeFiles/unit_tests.dir/flags.make
check/CMakeFiles/unit_tests.dir/TestIpx.cpp.o: ../check/TestIpx.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object check/CMakeFiles/unit_tests.dir/TestIpx.cpp.o"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/unit_tests.dir/TestIpx.cpp.o -c /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestIpx.cpp

check/CMakeFiles/unit_tests.dir/TestIpx.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/unit_tests.dir/TestIpx.cpp.i"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestIpx.cpp > CMakeFiles/unit_tests.dir/TestIpx.cpp.i

check/CMakeFiles/unit_tests.dir/TestIpx.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/unit_tests.dir/TestIpx.cpp.s"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/edeakins/Thesis/EQLPSolver/DHiGHS/check/TestIpx.cpp -o CMakeFiles/unit_tests.dir/TestIpx.cpp.s

check/CMakeFiles/unit_tests.dir/TestIpx.cpp.o.requires:

.PHONY : check/CMakeFiles/unit_tests.dir/TestIpx.cpp.o.requires

check/CMakeFiles/unit_tests.dir/TestIpx.cpp.o.provides: check/CMakeFiles/unit_tests.dir/TestIpx.cpp.o.requires
	$(MAKE) -f check/CMakeFiles/unit_tests.dir/build.make check/CMakeFiles/unit_tests.dir/TestIpx.cpp.o.provides.build
.PHONY : check/CMakeFiles/unit_tests.dir/TestIpx.cpp.o.provides

check/CMakeFiles/unit_tests.dir/TestIpx.cpp.o.provides.build: check/CMakeFiles/unit_tests.dir/TestIpx.cpp.o


# Object files for target unit_tests
unit_tests_OBJECTS = \
"CMakeFiles/unit_tests.dir/TestMain.cpp.o" \
"CMakeFiles/unit_tests.dir/TestOptions.cpp.o" \
"CMakeFiles/unit_tests.dir/TestPresolve.cpp.o" \
"CMakeFiles/unit_tests.dir/TestIO.cpp.o" \
"CMakeFiles/unit_tests.dir/TestSort.cpp.o" \
"CMakeFiles/unit_tests.dir/TestSetup.cpp.o" \
"CMakeFiles/unit_tests.dir/TestFilereader.cpp.o" \
"CMakeFiles/unit_tests.dir/TestInfo.cpp.o" \
"CMakeFiles/unit_tests.dir/TestBasisSolves.cpp.o" \
"CMakeFiles/unit_tests.dir/TestLpValidation.cpp.o" \
"CMakeFiles/unit_tests.dir/TestLpModification.cpp.o" \
"CMakeFiles/unit_tests.dir/Avgas.cpp.o" \
"CMakeFiles/unit_tests.dir/TestIpx.cpp.o"

# External object files for target unit_tests
unit_tests_EXTERNAL_OBJECTS =

bin/unit_tests: check/CMakeFiles/unit_tests.dir/TestMain.cpp.o
bin/unit_tests: check/CMakeFiles/unit_tests.dir/TestOptions.cpp.o
bin/unit_tests: check/CMakeFiles/unit_tests.dir/TestPresolve.cpp.o
bin/unit_tests: check/CMakeFiles/unit_tests.dir/TestIO.cpp.o
bin/unit_tests: check/CMakeFiles/unit_tests.dir/TestSort.cpp.o
bin/unit_tests: check/CMakeFiles/unit_tests.dir/TestSetup.cpp.o
bin/unit_tests: check/CMakeFiles/unit_tests.dir/TestFilereader.cpp.o
bin/unit_tests: check/CMakeFiles/unit_tests.dir/TestInfo.cpp.o
bin/unit_tests: check/CMakeFiles/unit_tests.dir/TestBasisSolves.cpp.o
bin/unit_tests: check/CMakeFiles/unit_tests.dir/TestLpValidation.cpp.o
bin/unit_tests: check/CMakeFiles/unit_tests.dir/TestLpModification.cpp.o
bin/unit_tests: check/CMakeFiles/unit_tests.dir/Avgas.cpp.o
bin/unit_tests: check/CMakeFiles/unit_tests.dir/TestIpx.cpp.o
bin/unit_tests: check/CMakeFiles/unit_tests.dir/build.make
bin/unit_tests: lib/libhighs.so.1.0.0
bin/unit_tests: lib/libipx.so
bin/unit_tests: lib/libbasiclu.so
bin/unit_tests: check/CMakeFiles/unit_tests.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Linking CXX executable ../bin/unit_tests"
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/unit_tests.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
check/CMakeFiles/unit_tests.dir/build: bin/unit_tests

.PHONY : check/CMakeFiles/unit_tests.dir/build

check/CMakeFiles/unit_tests.dir/requires: check/CMakeFiles/unit_tests.dir/TestMain.cpp.o.requires
check/CMakeFiles/unit_tests.dir/requires: check/CMakeFiles/unit_tests.dir/TestOptions.cpp.o.requires
check/CMakeFiles/unit_tests.dir/requires: check/CMakeFiles/unit_tests.dir/TestPresolve.cpp.o.requires
check/CMakeFiles/unit_tests.dir/requires: check/CMakeFiles/unit_tests.dir/TestIO.cpp.o.requires
check/CMakeFiles/unit_tests.dir/requires: check/CMakeFiles/unit_tests.dir/TestSort.cpp.o.requires
check/CMakeFiles/unit_tests.dir/requires: check/CMakeFiles/unit_tests.dir/TestSetup.cpp.o.requires
check/CMakeFiles/unit_tests.dir/requires: check/CMakeFiles/unit_tests.dir/TestFilereader.cpp.o.requires
check/CMakeFiles/unit_tests.dir/requires: check/CMakeFiles/unit_tests.dir/TestInfo.cpp.o.requires
check/CMakeFiles/unit_tests.dir/requires: check/CMakeFiles/unit_tests.dir/TestBasisSolves.cpp.o.requires
check/CMakeFiles/unit_tests.dir/requires: check/CMakeFiles/unit_tests.dir/TestLpValidation.cpp.o.requires
check/CMakeFiles/unit_tests.dir/requires: check/CMakeFiles/unit_tests.dir/TestLpModification.cpp.o.requires
check/CMakeFiles/unit_tests.dir/requires: check/CMakeFiles/unit_tests.dir/Avgas.cpp.o.requires
check/CMakeFiles/unit_tests.dir/requires: check/CMakeFiles/unit_tests.dir/TestIpx.cpp.o.requires

.PHONY : check/CMakeFiles/unit_tests.dir/requires

check/CMakeFiles/unit_tests.dir/clean:
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check && $(CMAKE_COMMAND) -P CMakeFiles/unit_tests.dir/cmake_clean.cmake
.PHONY : check/CMakeFiles/unit_tests.dir/clean

check/CMakeFiles/unit_tests.dir/depend:
	cd /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/edeakins/Thesis/EQLPSolver/DHiGHS /home/edeakins/Thesis/EQLPSolver/DHiGHS/check /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check /home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild/check/CMakeFiles/unit_tests.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : check/CMakeFiles/unit_tests.dir/depend


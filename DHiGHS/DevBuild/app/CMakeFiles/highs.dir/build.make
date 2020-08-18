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

# Include any dependencies generated for this target.
include app/CMakeFiles/highs.dir/depend.make

# Include the progress variables for this target.
include app/CMakeFiles/highs.dir/progress.make

# Include the compile flags for this target's objects.
include app/CMakeFiles/highs.dir/flags.make

app/CMakeFiles/highs.dir/RunHighs.cpp.o: app/CMakeFiles/highs.dir/flags.make
app/CMakeFiles/highs.dir/RunHighs.cpp.o: ../app/RunHighs.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/edeakins/EQLPSolver/DHiGHS/DevBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object app/CMakeFiles/highs.dir/RunHighs.cpp.o"
	cd /home/edeakins/EQLPSolver/DHiGHS/DevBuild/app && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/highs.dir/RunHighs.cpp.o -c /home/edeakins/EQLPSolver/DHiGHS/app/RunHighs.cpp

app/CMakeFiles/highs.dir/RunHighs.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/highs.dir/RunHighs.cpp.i"
	cd /home/edeakins/EQLPSolver/DHiGHS/DevBuild/app && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/edeakins/EQLPSolver/DHiGHS/app/RunHighs.cpp > CMakeFiles/highs.dir/RunHighs.cpp.i

app/CMakeFiles/highs.dir/RunHighs.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/highs.dir/RunHighs.cpp.s"
	cd /home/edeakins/EQLPSolver/DHiGHS/DevBuild/app && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/edeakins/EQLPSolver/DHiGHS/app/RunHighs.cpp -o CMakeFiles/highs.dir/RunHighs.cpp.s

# Object files for target highs
highs_OBJECTS = \
"CMakeFiles/highs.dir/RunHighs.cpp.o"

# External object files for target highs
highs_EXTERNAL_OBJECTS =

bin/highs: app/CMakeFiles/highs.dir/RunHighs.cpp.o
bin/highs: app/CMakeFiles/highs.dir/build.make
bin/highs: lib/libhighs.so.1.0.0
bin/highs: lib/libipx.so
bin/highs: lib/libbasiclu.so
bin/highs: app/CMakeFiles/highs.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/edeakins/EQLPSolver/DHiGHS/DevBuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/highs"
	cd /home/edeakins/EQLPSolver/DHiGHS/DevBuild/app && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/highs.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
app/CMakeFiles/highs.dir/build: bin/highs

.PHONY : app/CMakeFiles/highs.dir/build

app/CMakeFiles/highs.dir/clean:
	cd /home/edeakins/EQLPSolver/DHiGHS/DevBuild/app && $(CMAKE_COMMAND) -P CMakeFiles/highs.dir/cmake_clean.cmake
.PHONY : app/CMakeFiles/highs.dir/clean

app/CMakeFiles/highs.dir/depend:
	cd /home/edeakins/EQLPSolver/DHiGHS/DevBuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/edeakins/EQLPSolver/DHiGHS /home/edeakins/EQLPSolver/DHiGHS/app /home/edeakins/EQLPSolver/DHiGHS/DevBuild /home/edeakins/EQLPSolver/DHiGHS/DevBuild/app /home/edeakins/EQLPSolver/DHiGHS/DevBuild/app/CMakeFiles/highs.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : app/CMakeFiles/highs.dir/depend


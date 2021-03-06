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
CMAKE_COMMAND = /home/maytee/software/clion-2018.1/bin/cmake/bin/cmake

# The command to remove a file.
RM = /home/maytee/software/clion-2018.1/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/maytee/Documents/2.29/finiteVolumeSolver

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug

# Include any dependencies generated for this target.
include src/eigen/bench/spbench/CMakeFiles/spsolver.dir/depend.make

# Include the progress variables for this target.
include src/eigen/bench/spbench/CMakeFiles/spsolver.dir/progress.make

# Include the compile flags for this target's objects.
include src/eigen/bench/spbench/CMakeFiles/spsolver.dir/flags.make

src/eigen/bench/spbench/CMakeFiles/spsolver.dir/sp_solver.cpp.o: src/eigen/bench/spbench/CMakeFiles/spsolver.dir/flags.make
src/eigen/bench/spbench/CMakeFiles/spsolver.dir/sp_solver.cpp.o: ../src/eigen/bench/spbench/sp_solver.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/eigen/bench/spbench/CMakeFiles/spsolver.dir/sp_solver.cpp.o"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/eigen/bench/spbench && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/spsolver.dir/sp_solver.cpp.o -c /home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/bench/spbench/sp_solver.cpp

src/eigen/bench/spbench/CMakeFiles/spsolver.dir/sp_solver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/spsolver.dir/sp_solver.cpp.i"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/eigen/bench/spbench && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/bench/spbench/sp_solver.cpp > CMakeFiles/spsolver.dir/sp_solver.cpp.i

src/eigen/bench/spbench/CMakeFiles/spsolver.dir/sp_solver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/spsolver.dir/sp_solver.cpp.s"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/eigen/bench/spbench && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/bench/spbench/sp_solver.cpp -o CMakeFiles/spsolver.dir/sp_solver.cpp.s

src/eigen/bench/spbench/CMakeFiles/spsolver.dir/sp_solver.cpp.o.requires:

.PHONY : src/eigen/bench/spbench/CMakeFiles/spsolver.dir/sp_solver.cpp.o.requires

src/eigen/bench/spbench/CMakeFiles/spsolver.dir/sp_solver.cpp.o.provides: src/eigen/bench/spbench/CMakeFiles/spsolver.dir/sp_solver.cpp.o.requires
	$(MAKE) -f src/eigen/bench/spbench/CMakeFiles/spsolver.dir/build.make src/eigen/bench/spbench/CMakeFiles/spsolver.dir/sp_solver.cpp.o.provides.build
.PHONY : src/eigen/bench/spbench/CMakeFiles/spsolver.dir/sp_solver.cpp.o.provides

src/eigen/bench/spbench/CMakeFiles/spsolver.dir/sp_solver.cpp.o.provides.build: src/eigen/bench/spbench/CMakeFiles/spsolver.dir/sp_solver.cpp.o


# Object files for target spsolver
spsolver_OBJECTS = \
"CMakeFiles/spsolver.dir/sp_solver.cpp.o"

# External object files for target spsolver
spsolver_EXTERNAL_OBJECTS =

src/eigen/bench/spbench/spsolver: src/eigen/bench/spbench/CMakeFiles/spsolver.dir/sp_solver.cpp.o
src/eigen/bench/spbench/spsolver: src/eigen/bench/spbench/CMakeFiles/spsolver.dir/build.make
src/eigen/bench/spbench/spsolver: /usr/lib/x86_64-linux-gnu/libsuperlu.so
src/eigen/bench/spbench/spsolver: src/eigen/blas/libeigen_blas_static.a
src/eigen/bench/spbench/spsolver: /usr/lib/x86_64-linux-gnu/librt.so
src/eigen/bench/spbench/spsolver: src/eigen/bench/spbench/CMakeFiles/spsolver.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable spsolver"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/eigen/bench/spbench && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/spsolver.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/eigen/bench/spbench/CMakeFiles/spsolver.dir/build: src/eigen/bench/spbench/spsolver

.PHONY : src/eigen/bench/spbench/CMakeFiles/spsolver.dir/build

src/eigen/bench/spbench/CMakeFiles/spsolver.dir/requires: src/eigen/bench/spbench/CMakeFiles/spsolver.dir/sp_solver.cpp.o.requires

.PHONY : src/eigen/bench/spbench/CMakeFiles/spsolver.dir/requires

src/eigen/bench/spbench/CMakeFiles/spsolver.dir/clean:
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/eigen/bench/spbench && $(CMAKE_COMMAND) -P CMakeFiles/spsolver.dir/cmake_clean.cmake
.PHONY : src/eigen/bench/spbench/CMakeFiles/spsolver.dir/clean

src/eigen/bench/spbench/CMakeFiles/spsolver.dir/depend:
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/maytee/Documents/2.29/finiteVolumeSolver /home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/bench/spbench /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/eigen/bench/spbench /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/eigen/bench/spbench/CMakeFiles/spsolver.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/eigen/bench/spbench/CMakeFiles/spsolver.dir/depend


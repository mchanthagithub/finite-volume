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
include src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/depend.make

# Include the progress variables for this target.
include src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/progress.make

# Include the compile flags for this target's objects.
include src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/flags.make

src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/dblat3.f.o: src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/flags.make
src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/dblat3.f.o: ../src/include/eigen/blas/testing/dblat3.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/dblat3.f.o"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/include/eigen/blas/testing && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/maytee/Documents/2.29/finiteVolumeSolver/src/include/eigen/blas/testing/dblat3.f -o CMakeFiles/dblat3.dir/dblat3.f.o

src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/dblat3.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/dblat3.dir/dblat3.f.i"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/include/eigen/blas/testing && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/maytee/Documents/2.29/finiteVolumeSolver/src/include/eigen/blas/testing/dblat3.f > CMakeFiles/dblat3.dir/dblat3.f.i

src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/dblat3.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/dblat3.dir/dblat3.f.s"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/include/eigen/blas/testing && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/maytee/Documents/2.29/finiteVolumeSolver/src/include/eigen/blas/testing/dblat3.f -o CMakeFiles/dblat3.dir/dblat3.f.s

src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/dblat3.f.o.requires:

.PHONY : src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/dblat3.f.o.requires

src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/dblat3.f.o.provides: src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/dblat3.f.o.requires
	$(MAKE) -f src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/build.make src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/dblat3.f.o.provides.build
.PHONY : src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/dblat3.f.o.provides

src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/dblat3.f.o.provides.build: src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/dblat3.f.o


# Object files for target dblat3
dblat3_OBJECTS = \
"CMakeFiles/dblat3.dir/dblat3.f.o"

# External object files for target dblat3
dblat3_EXTERNAL_OBJECTS =

src/include/eigen/blas/testing/dblat3: src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/dblat3.f.o
src/include/eigen/blas/testing/dblat3: src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/build.make
src/include/eigen/blas/testing/dblat3: src/include/eigen/blas/libeigen_blas.so
src/include/eigen/blas/testing/dblat3: src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable dblat3"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/include/eigen/blas/testing && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dblat3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/build: src/include/eigen/blas/testing/dblat3

.PHONY : src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/build

src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/requires: src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/dblat3.f.o.requires

.PHONY : src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/requires

src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/clean:
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/include/eigen/blas/testing && $(CMAKE_COMMAND) -P CMakeFiles/dblat3.dir/cmake_clean.cmake
.PHONY : src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/clean

src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/depend:
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/maytee/Documents/2.29/finiteVolumeSolver /home/maytee/Documents/2.29/finiteVolumeSolver/src/include/eigen/blas/testing /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/include/eigen/blas/testing /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/include/eigen/blas/testing/CMakeFiles/dblat3.dir/depend


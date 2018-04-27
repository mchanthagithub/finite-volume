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
include src/eigen/blas/testing/CMakeFiles/sblat3.dir/depend.make

# Include the progress variables for this target.
include src/eigen/blas/testing/CMakeFiles/sblat3.dir/progress.make

# Include the compile flags for this target's objects.
include src/eigen/blas/testing/CMakeFiles/sblat3.dir/flags.make

src/eigen/blas/testing/CMakeFiles/sblat3.dir/sblat3.f.o: src/eigen/blas/testing/CMakeFiles/sblat3.dir/flags.make
src/eigen/blas/testing/CMakeFiles/sblat3.dir/sblat3.f.o: ../src/eigen/blas/testing/sblat3.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object src/eigen/blas/testing/CMakeFiles/sblat3.dir/sblat3.f.o"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/eigen/blas/testing && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/blas/testing/sblat3.f -o CMakeFiles/sblat3.dir/sblat3.f.o

src/eigen/blas/testing/CMakeFiles/sblat3.dir/sblat3.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/sblat3.dir/sblat3.f.i"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/eigen/blas/testing && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/blas/testing/sblat3.f > CMakeFiles/sblat3.dir/sblat3.f.i

src/eigen/blas/testing/CMakeFiles/sblat3.dir/sblat3.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/sblat3.dir/sblat3.f.s"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/eigen/blas/testing && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/blas/testing/sblat3.f -o CMakeFiles/sblat3.dir/sblat3.f.s

src/eigen/blas/testing/CMakeFiles/sblat3.dir/sblat3.f.o.requires:

.PHONY : src/eigen/blas/testing/CMakeFiles/sblat3.dir/sblat3.f.o.requires

src/eigen/blas/testing/CMakeFiles/sblat3.dir/sblat3.f.o.provides: src/eigen/blas/testing/CMakeFiles/sblat3.dir/sblat3.f.o.requires
	$(MAKE) -f src/eigen/blas/testing/CMakeFiles/sblat3.dir/build.make src/eigen/blas/testing/CMakeFiles/sblat3.dir/sblat3.f.o.provides.build
.PHONY : src/eigen/blas/testing/CMakeFiles/sblat3.dir/sblat3.f.o.provides

src/eigen/blas/testing/CMakeFiles/sblat3.dir/sblat3.f.o.provides.build: src/eigen/blas/testing/CMakeFiles/sblat3.dir/sblat3.f.o


# Object files for target sblat3
sblat3_OBJECTS = \
"CMakeFiles/sblat3.dir/sblat3.f.o"

# External object files for target sblat3
sblat3_EXTERNAL_OBJECTS =

src/eigen/blas/testing/sblat3: src/eigen/blas/testing/CMakeFiles/sblat3.dir/sblat3.f.o
src/eigen/blas/testing/sblat3: src/eigen/blas/testing/CMakeFiles/sblat3.dir/build.make
src/eigen/blas/testing/sblat3: src/eigen/blas/libeigen_blas.so
src/eigen/blas/testing/sblat3: src/eigen/blas/testing/CMakeFiles/sblat3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable sblat3"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/eigen/blas/testing && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sblat3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/eigen/blas/testing/CMakeFiles/sblat3.dir/build: src/eigen/blas/testing/sblat3

.PHONY : src/eigen/blas/testing/CMakeFiles/sblat3.dir/build

src/eigen/blas/testing/CMakeFiles/sblat3.dir/requires: src/eigen/blas/testing/CMakeFiles/sblat3.dir/sblat3.f.o.requires

.PHONY : src/eigen/blas/testing/CMakeFiles/sblat3.dir/requires

src/eigen/blas/testing/CMakeFiles/sblat3.dir/clean:
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/eigen/blas/testing && $(CMAKE_COMMAND) -P CMakeFiles/sblat3.dir/cmake_clean.cmake
.PHONY : src/eigen/blas/testing/CMakeFiles/sblat3.dir/clean

src/eigen/blas/testing/CMakeFiles/sblat3.dir/depend:
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/maytee/Documents/2.29/finiteVolumeSolver /home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/blas/testing /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/eigen/blas/testing /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/eigen/blas/testing/CMakeFiles/sblat3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/eigen/blas/testing/CMakeFiles/sblat3.dir/depend


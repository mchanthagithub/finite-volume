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
include src/CMakeFiles/finiteVolumeSolver.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/finiteVolumeSolver.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/finiteVolumeSolver.dir/flags.make

src/CMakeFiles/finiteVolumeSolver.dir/main.cpp.o: src/CMakeFiles/finiteVolumeSolver.dir/flags.make
src/CMakeFiles/finiteVolumeSolver.dir/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/finiteVolumeSolver.dir/main.cpp.o"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/finiteVolumeSolver.dir/main.cpp.o -c /home/maytee/Documents/2.29/finiteVolumeSolver/src/main.cpp

src/CMakeFiles/finiteVolumeSolver.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/finiteVolumeSolver.dir/main.cpp.i"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/maytee/Documents/2.29/finiteVolumeSolver/src/main.cpp > CMakeFiles/finiteVolumeSolver.dir/main.cpp.i

src/CMakeFiles/finiteVolumeSolver.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/finiteVolumeSolver.dir/main.cpp.s"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/maytee/Documents/2.29/finiteVolumeSolver/src/main.cpp -o CMakeFiles/finiteVolumeSolver.dir/main.cpp.s

src/CMakeFiles/finiteVolumeSolver.dir/main.cpp.o.requires:

.PHONY : src/CMakeFiles/finiteVolumeSolver.dir/main.cpp.o.requires

src/CMakeFiles/finiteVolumeSolver.dir/main.cpp.o.provides: src/CMakeFiles/finiteVolumeSolver.dir/main.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/finiteVolumeSolver.dir/build.make src/CMakeFiles/finiteVolumeSolver.dir/main.cpp.o.provides.build
.PHONY : src/CMakeFiles/finiteVolumeSolver.dir/main.cpp.o.provides

src/CMakeFiles/finiteVolumeSolver.dir/main.cpp.o.provides.build: src/CMakeFiles/finiteVolumeSolver.dir/main.cpp.o


src/CMakeFiles/finiteVolumeSolver.dir/Grids/Grid.cpp.o: src/CMakeFiles/finiteVolumeSolver.dir/flags.make
src/CMakeFiles/finiteVolumeSolver.dir/Grids/Grid.cpp.o: ../src/Grids/Grid.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/CMakeFiles/finiteVolumeSolver.dir/Grids/Grid.cpp.o"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/finiteVolumeSolver.dir/Grids/Grid.cpp.o -c /home/maytee/Documents/2.29/finiteVolumeSolver/src/Grids/Grid.cpp

src/CMakeFiles/finiteVolumeSolver.dir/Grids/Grid.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/finiteVolumeSolver.dir/Grids/Grid.cpp.i"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/maytee/Documents/2.29/finiteVolumeSolver/src/Grids/Grid.cpp > CMakeFiles/finiteVolumeSolver.dir/Grids/Grid.cpp.i

src/CMakeFiles/finiteVolumeSolver.dir/Grids/Grid.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/finiteVolumeSolver.dir/Grids/Grid.cpp.s"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/maytee/Documents/2.29/finiteVolumeSolver/src/Grids/Grid.cpp -o CMakeFiles/finiteVolumeSolver.dir/Grids/Grid.cpp.s

src/CMakeFiles/finiteVolumeSolver.dir/Grids/Grid.cpp.o.requires:

.PHONY : src/CMakeFiles/finiteVolumeSolver.dir/Grids/Grid.cpp.o.requires

src/CMakeFiles/finiteVolumeSolver.dir/Grids/Grid.cpp.o.provides: src/CMakeFiles/finiteVolumeSolver.dir/Grids/Grid.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/finiteVolumeSolver.dir/build.make src/CMakeFiles/finiteVolumeSolver.dir/Grids/Grid.cpp.o.provides.build
.PHONY : src/CMakeFiles/finiteVolumeSolver.dir/Grids/Grid.cpp.o.provides

src/CMakeFiles/finiteVolumeSolver.dir/Grids/Grid.cpp.o.provides.build: src/CMakeFiles/finiteVolumeSolver.dir/Grids/Grid.cpp.o


src/CMakeFiles/finiteVolumeSolver.dir/Grids/CartesianGrid.cpp.o: src/CMakeFiles/finiteVolumeSolver.dir/flags.make
src/CMakeFiles/finiteVolumeSolver.dir/Grids/CartesianGrid.cpp.o: ../src/Grids/CartesianGrid.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/CMakeFiles/finiteVolumeSolver.dir/Grids/CartesianGrid.cpp.o"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/finiteVolumeSolver.dir/Grids/CartesianGrid.cpp.o -c /home/maytee/Documents/2.29/finiteVolumeSolver/src/Grids/CartesianGrid.cpp

src/CMakeFiles/finiteVolumeSolver.dir/Grids/CartesianGrid.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/finiteVolumeSolver.dir/Grids/CartesianGrid.cpp.i"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/maytee/Documents/2.29/finiteVolumeSolver/src/Grids/CartesianGrid.cpp > CMakeFiles/finiteVolumeSolver.dir/Grids/CartesianGrid.cpp.i

src/CMakeFiles/finiteVolumeSolver.dir/Grids/CartesianGrid.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/finiteVolumeSolver.dir/Grids/CartesianGrid.cpp.s"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/maytee/Documents/2.29/finiteVolumeSolver/src/Grids/CartesianGrid.cpp -o CMakeFiles/finiteVolumeSolver.dir/Grids/CartesianGrid.cpp.s

src/CMakeFiles/finiteVolumeSolver.dir/Grids/CartesianGrid.cpp.o.requires:

.PHONY : src/CMakeFiles/finiteVolumeSolver.dir/Grids/CartesianGrid.cpp.o.requires

src/CMakeFiles/finiteVolumeSolver.dir/Grids/CartesianGrid.cpp.o.provides: src/CMakeFiles/finiteVolumeSolver.dir/Grids/CartesianGrid.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/finiteVolumeSolver.dir/build.make src/CMakeFiles/finiteVolumeSolver.dir/Grids/CartesianGrid.cpp.o.provides.build
.PHONY : src/CMakeFiles/finiteVolumeSolver.dir/Grids/CartesianGrid.cpp.o.provides

src/CMakeFiles/finiteVolumeSolver.dir/Grids/CartesianGrid.cpp.o.provides.build: src/CMakeFiles/finiteVolumeSolver.dir/Grids/CartesianGrid.cpp.o


# Object files for target finiteVolumeSolver
finiteVolumeSolver_OBJECTS = \
"CMakeFiles/finiteVolumeSolver.dir/main.cpp.o" \
"CMakeFiles/finiteVolumeSolver.dir/Grids/Grid.cpp.o" \
"CMakeFiles/finiteVolumeSolver.dir/Grids/CartesianGrid.cpp.o"

# External object files for target finiteVolumeSolver
finiteVolumeSolver_EXTERNAL_OBJECTS =

bin/finiteVolumeSolver: src/CMakeFiles/finiteVolumeSolver.dir/main.cpp.o
bin/finiteVolumeSolver: src/CMakeFiles/finiteVolumeSolver.dir/Grids/Grid.cpp.o
bin/finiteVolumeSolver: src/CMakeFiles/finiteVolumeSolver.dir/Grids/CartesianGrid.cpp.o
bin/finiteVolumeSolver: src/CMakeFiles/finiteVolumeSolver.dir/build.make
bin/finiteVolumeSolver: src/CMakeFiles/finiteVolumeSolver.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable ../bin/finiteVolumeSolver"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/finiteVolumeSolver.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/finiteVolumeSolver.dir/build: bin/finiteVolumeSolver

.PHONY : src/CMakeFiles/finiteVolumeSolver.dir/build

src/CMakeFiles/finiteVolumeSolver.dir/requires: src/CMakeFiles/finiteVolumeSolver.dir/main.cpp.o.requires
src/CMakeFiles/finiteVolumeSolver.dir/requires: src/CMakeFiles/finiteVolumeSolver.dir/Grids/Grid.cpp.o.requires
src/CMakeFiles/finiteVolumeSolver.dir/requires: src/CMakeFiles/finiteVolumeSolver.dir/Grids/CartesianGrid.cpp.o.requires

.PHONY : src/CMakeFiles/finiteVolumeSolver.dir/requires

src/CMakeFiles/finiteVolumeSolver.dir/clean:
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src && $(CMAKE_COMMAND) -P CMakeFiles/finiteVolumeSolver.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/finiteVolumeSolver.dir/clean

src/CMakeFiles/finiteVolumeSolver.dir/depend:
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/maytee/Documents/2.29/finiteVolumeSolver /home/maytee/Documents/2.29/finiteVolumeSolver/src /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/CMakeFiles/finiteVolumeSolver.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/finiteVolumeSolver.dir/depend


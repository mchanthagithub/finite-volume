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
include src/eigen/test/CMakeFiles/sparse_block_4.dir/depend.make

# Include the progress variables for this target.
include src/eigen/test/CMakeFiles/sparse_block_4.dir/progress.make

# Include the compile flags for this target's objects.
include src/eigen/test/CMakeFiles/sparse_block_4.dir/flags.make

src/eigen/test/CMakeFiles/sparse_block_4.dir/sparse_block.cpp.o: src/eigen/test/CMakeFiles/sparse_block_4.dir/flags.make
src/eigen/test/CMakeFiles/sparse_block_4.dir/sparse_block.cpp.o: ../src/eigen/test/sparse_block.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/eigen/test/CMakeFiles/sparse_block_4.dir/sparse_block.cpp.o"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/eigen/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/sparse_block_4.dir/sparse_block.cpp.o -c /home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/test/sparse_block.cpp

src/eigen/test/CMakeFiles/sparse_block_4.dir/sparse_block.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sparse_block_4.dir/sparse_block.cpp.i"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/eigen/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/test/sparse_block.cpp > CMakeFiles/sparse_block_4.dir/sparse_block.cpp.i

src/eigen/test/CMakeFiles/sparse_block_4.dir/sparse_block.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sparse_block_4.dir/sparse_block.cpp.s"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/eigen/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/test/sparse_block.cpp -o CMakeFiles/sparse_block_4.dir/sparse_block.cpp.s

src/eigen/test/CMakeFiles/sparse_block_4.dir/sparse_block.cpp.o.requires:

.PHONY : src/eigen/test/CMakeFiles/sparse_block_4.dir/sparse_block.cpp.o.requires

src/eigen/test/CMakeFiles/sparse_block_4.dir/sparse_block.cpp.o.provides: src/eigen/test/CMakeFiles/sparse_block_4.dir/sparse_block.cpp.o.requires
	$(MAKE) -f src/eigen/test/CMakeFiles/sparse_block_4.dir/build.make src/eigen/test/CMakeFiles/sparse_block_4.dir/sparse_block.cpp.o.provides.build
.PHONY : src/eigen/test/CMakeFiles/sparse_block_4.dir/sparse_block.cpp.o.provides

src/eigen/test/CMakeFiles/sparse_block_4.dir/sparse_block.cpp.o.provides.build: src/eigen/test/CMakeFiles/sparse_block_4.dir/sparse_block.cpp.o


# Object files for target sparse_block_4
sparse_block_4_OBJECTS = \
"CMakeFiles/sparse_block_4.dir/sparse_block.cpp.o"

# External object files for target sparse_block_4
sparse_block_4_EXTERNAL_OBJECTS =

src/eigen/test/sparse_block_4: src/eigen/test/CMakeFiles/sparse_block_4.dir/sparse_block.cpp.o
src/eigen/test/sparse_block_4: src/eigen/test/CMakeFiles/sparse_block_4.dir/build.make
src/eigen/test/sparse_block_4: src/eigen/test/CMakeFiles/sparse_block_4.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable sparse_block_4"
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/eigen/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sparse_block_4.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/eigen/test/CMakeFiles/sparse_block_4.dir/build: src/eigen/test/sparse_block_4

.PHONY : src/eigen/test/CMakeFiles/sparse_block_4.dir/build

src/eigen/test/CMakeFiles/sparse_block_4.dir/requires: src/eigen/test/CMakeFiles/sparse_block_4.dir/sparse_block.cpp.o.requires

.PHONY : src/eigen/test/CMakeFiles/sparse_block_4.dir/requires

src/eigen/test/CMakeFiles/sparse_block_4.dir/clean:
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/eigen/test && $(CMAKE_COMMAND) -P CMakeFiles/sparse_block_4.dir/cmake_clean.cmake
.PHONY : src/eigen/test/CMakeFiles/sparse_block_4.dir/clean

src/eigen/test/CMakeFiles/sparse_block_4.dir/depend:
	cd /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/maytee/Documents/2.29/finiteVolumeSolver /home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/test /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/eigen/test /home/maytee/Documents/2.29/finiteVolumeSolver/cmake-build-debug/src/eigen/test/CMakeFiles/sparse_block_4.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/eigen/test/CMakeFiles/sparse_block_4.dir/depend


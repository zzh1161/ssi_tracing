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
CMAKE_SOURCE_DIR = /home/zhang-zh/workspace/CAD/SSI_tracing

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/zhang-zh/workspace/CAD/SSI_tracing/build

# Include any dependencies generated for this target.
include test/CMakeFiles/test_tn.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/test_tn.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/test_tn.dir/flags.make

test/CMakeFiles/test_tn.dir/test_tn.cpp.o: test/CMakeFiles/test_tn.dir/flags.make
test/CMakeFiles/test_tn.dir/test_tn.cpp.o: ../test/test_tn.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zhang-zh/workspace/CAD/SSI_tracing/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/test_tn.dir/test_tn.cpp.o"
	cd /home/zhang-zh/workspace/CAD/SSI_tracing/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_tn.dir/test_tn.cpp.o -c /home/zhang-zh/workspace/CAD/SSI_tracing/test/test_tn.cpp

test/CMakeFiles/test_tn.dir/test_tn.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_tn.dir/test_tn.cpp.i"
	cd /home/zhang-zh/workspace/CAD/SSI_tracing/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zhang-zh/workspace/CAD/SSI_tracing/test/test_tn.cpp > CMakeFiles/test_tn.dir/test_tn.cpp.i

test/CMakeFiles/test_tn.dir/test_tn.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_tn.dir/test_tn.cpp.s"
	cd /home/zhang-zh/workspace/CAD/SSI_tracing/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zhang-zh/workspace/CAD/SSI_tracing/test/test_tn.cpp -o CMakeFiles/test_tn.dir/test_tn.cpp.s

# Object files for target test_tn
test_tn_OBJECTS = \
"CMakeFiles/test_tn.dir/test_tn.cpp.o"

# External object files for target test_tn
test_tn_EXTERNAL_OBJECTS =

test/test_tn: test/CMakeFiles/test_tn.dir/test_tn.cpp.o
test/test_tn: test/CMakeFiles/test_tn.dir/build.make
test/test_tn: test/CMakeFiles/test_tn.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/zhang-zh/workspace/CAD/SSI_tracing/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test_tn"
	cd /home/zhang-zh/workspace/CAD/SSI_tracing/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_tn.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/test_tn.dir/build: test/test_tn

.PHONY : test/CMakeFiles/test_tn.dir/build

test/CMakeFiles/test_tn.dir/clean:
	cd /home/zhang-zh/workspace/CAD/SSI_tracing/build/test && $(CMAKE_COMMAND) -P CMakeFiles/test_tn.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/test_tn.dir/clean

test/CMakeFiles/test_tn.dir/depend:
	cd /home/zhang-zh/workspace/CAD/SSI_tracing/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zhang-zh/workspace/CAD/SSI_tracing /home/zhang-zh/workspace/CAD/SSI_tracing/test /home/zhang-zh/workspace/CAD/SSI_tracing/build /home/zhang-zh/workspace/CAD/SSI_tracing/build/test /home/zhang-zh/workspace/CAD/SSI_tracing/build/test/CMakeFiles/test_tn.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/test_tn.dir/depend

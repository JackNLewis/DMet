# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/jacklewis/Documents/work/year3/DMet

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/jacklewis/Documents/work/year3/DMet/build

# Include any dependencies generated for this target.
include test/CMakeFiles/Google_tests_run.dir/depend.make
# Include the progress variables for this target.
include test/CMakeFiles/Google_tests_run.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/Google_tests_run.dir/flags.make

test/CMakeFiles/Google_tests_run.dir/DMet/PointDistanceTests.cpp.o: test/CMakeFiles/Google_tests_run.dir/flags.make
test/CMakeFiles/Google_tests_run.dir/DMet/PointDistanceTests.cpp.o: ../test/DMet/PointDistanceTests.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/jacklewis/Documents/work/year3/DMet/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/Google_tests_run.dir/DMet/PointDistanceTests.cpp.o"
	cd /Users/jacklewis/Documents/work/year3/DMet/build/test && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Google_tests_run.dir/DMet/PointDistanceTests.cpp.o -c /Users/jacklewis/Documents/work/year3/DMet/test/DMet/PointDistanceTests.cpp

test/CMakeFiles/Google_tests_run.dir/DMet/PointDistanceTests.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Google_tests_run.dir/DMet/PointDistanceTests.cpp.i"
	cd /Users/jacklewis/Documents/work/year3/DMet/build/test && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/jacklewis/Documents/work/year3/DMet/test/DMet/PointDistanceTests.cpp > CMakeFiles/Google_tests_run.dir/DMet/PointDistanceTests.cpp.i

test/CMakeFiles/Google_tests_run.dir/DMet/PointDistanceTests.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Google_tests_run.dir/DMet/PointDistanceTests.cpp.s"
	cd /Users/jacklewis/Documents/work/year3/DMet/build/test && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/jacklewis/Documents/work/year3/DMet/test/DMet/PointDistanceTests.cpp -o CMakeFiles/Google_tests_run.dir/DMet/PointDistanceTests.cpp.s

# Object files for target Google_tests_run
Google_tests_run_OBJECTS = \
"CMakeFiles/Google_tests_run.dir/DMet/PointDistanceTests.cpp.o"

# External object files for target Google_tests_run
Google_tests_run_EXTERNAL_OBJECTS =

test/Google_tests_run: test/CMakeFiles/Google_tests_run.dir/DMet/PointDistanceTests.cpp.o
test/Google_tests_run: test/CMakeFiles/Google_tests_run.dir/build.make
test/Google_tests_run: libDMet.dylib
test/Google_tests_run: lib/libgtest_maind.a
test/Google_tests_run: /usr/local/Cellar/gmp/6.2.1_1/lib/libgmp.a
test/Google_tests_run: /usr/local/Cellar/mpfr/4.1.0/lib/libmpfr.a
test/Google_tests_run: lib/libgtestd.a
test/Google_tests_run: test/CMakeFiles/Google_tests_run.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/jacklewis/Documents/work/year3/DMet/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Google_tests_run"
	cd /Users/jacklewis/Documents/work/year3/DMet/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Google_tests_run.dir/link.txt --verbose=$(VERBOSE)
	cd /Users/jacklewis/Documents/work/year3/DMet/build/test && /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -D TEST_TARGET=Google_tests_run -D TEST_EXECUTABLE=/Users/jacklewis/Documents/work/year3/DMet/build/test/Google_tests_run -D TEST_EXECUTOR= -D TEST_WORKING_DIR=/Users/jacklewis/Documents/work/year3/DMet/build/test -D TEST_EXTRA_ARGS= -D TEST_PROPERTIES= -D TEST_PREFIX= -D TEST_SUFFIX= -D NO_PRETTY_TYPES=FALSE -D NO_PRETTY_VALUES=FALSE -D TEST_LIST=Google_tests_run_TESTS -D CTEST_FILE=/Users/jacklewis/Documents/work/year3/DMet/build/test/Google_tests_run[1]_tests.cmake -D TEST_DISCOVERY_TIMEOUT=5 -D TEST_XML_OUTPUT_DIR= -P /Applications/CLion.app/Contents/bin/cmake/mac/share/cmake-3.20/Modules/GoogleTestAddTests.cmake

# Rule to build all files generated by this target.
test/CMakeFiles/Google_tests_run.dir/build: test/Google_tests_run
.PHONY : test/CMakeFiles/Google_tests_run.dir/build

test/CMakeFiles/Google_tests_run.dir/clean:
	cd /Users/jacklewis/Documents/work/year3/DMet/build/test && $(CMAKE_COMMAND) -P CMakeFiles/Google_tests_run.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/Google_tests_run.dir/clean

test/CMakeFiles/Google_tests_run.dir/depend:
	cd /Users/jacklewis/Documents/work/year3/DMet/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/jacklewis/Documents/work/year3/DMet /Users/jacklewis/Documents/work/year3/DMet/test /Users/jacklewis/Documents/work/year3/DMet/build /Users/jacklewis/Documents/work/year3/DMet/build/test /Users/jacklewis/Documents/work/year3/DMet/build/test/CMakeFiles/Google_tests_run.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/Google_tests_run.dir/depend


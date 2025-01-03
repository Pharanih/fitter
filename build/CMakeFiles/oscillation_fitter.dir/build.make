# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.30

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jake/Projects/Fitter/StatOnly

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jake/Projects/Fitter/StatOnly/build

# Include any dependencies generated for this target.
include CMakeFiles/oscillation_fitter.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/oscillation_fitter.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/oscillation_fitter.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/oscillation_fitter.dir/flags.make

CMakeFiles/oscillation_fitter.dir/oscillation_fitter.cxx.o: CMakeFiles/oscillation_fitter.dir/flags.make
CMakeFiles/oscillation_fitter.dir/oscillation_fitter.cxx.o: /home/jake/Projects/Fitter/StatOnly/oscillation_fitter.cxx
CMakeFiles/oscillation_fitter.dir/oscillation_fitter.cxx.o: CMakeFiles/oscillation_fitter.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/jake/Projects/Fitter/StatOnly/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/oscillation_fitter.dir/oscillation_fitter.cxx.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/oscillation_fitter.dir/oscillation_fitter.cxx.o -MF CMakeFiles/oscillation_fitter.dir/oscillation_fitter.cxx.o.d -o CMakeFiles/oscillation_fitter.dir/oscillation_fitter.cxx.o -c /home/jake/Projects/Fitter/StatOnly/oscillation_fitter.cxx

CMakeFiles/oscillation_fitter.dir/oscillation_fitter.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/oscillation_fitter.dir/oscillation_fitter.cxx.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jake/Projects/Fitter/StatOnly/oscillation_fitter.cxx > CMakeFiles/oscillation_fitter.dir/oscillation_fitter.cxx.i

CMakeFiles/oscillation_fitter.dir/oscillation_fitter.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/oscillation_fitter.dir/oscillation_fitter.cxx.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jake/Projects/Fitter/StatOnly/oscillation_fitter.cxx -o CMakeFiles/oscillation_fitter.dir/oscillation_fitter.cxx.s

# Object files for target oscillation_fitter
oscillation_fitter_OBJECTS = \
"CMakeFiles/oscillation_fitter.dir/oscillation_fitter.cxx.o"

# External object files for target oscillation_fitter
oscillation_fitter_EXTERNAL_OBJECTS =

oscillation_fitter: CMakeFiles/oscillation_fitter.dir/oscillation_fitter.cxx.o
oscillation_fitter: CMakeFiles/oscillation_fitter.dir/build.make
oscillation_fitter: /usr/lib64/root/libCore.so
oscillation_fitter: /usr/lib64/root/libImt.so
oscillation_fitter: /usr/lib64/root/libRIO.so
oscillation_fitter: /usr/lib64/root/libNet.so
oscillation_fitter: /usr/lib64/root/libHist.so
oscillation_fitter: /usr/lib64/root/libGraf.so
oscillation_fitter: /usr/lib64/root/libGraf3d.so
oscillation_fitter: /usr/lib64/root/libGpad.so
oscillation_fitter: /usr/lib64/root/libROOTDataFrame.so
oscillation_fitter: /usr/lib64/root/libTree.so
oscillation_fitter: /usr/lib64/root/libTreePlayer.so
oscillation_fitter: /usr/lib64/root/libRint.so
oscillation_fitter: /usr/lib64/root/libPostscript.so
oscillation_fitter: /usr/lib64/root/libMatrix.so
oscillation_fitter: /usr/lib64/root/libPhysics.so
oscillation_fitter: /usr/lib64/root/libMathCore.so
oscillation_fitter: /usr/lib64/root/libThread.so
oscillation_fitter: /usr/lib64/root/libMultiProc.so
oscillation_fitter: /usr/lib64/root/libROOTVecOps.so
oscillation_fitter: /usr/lib64/root/libMinuit.so
oscillation_fitter: CMakeFiles/oscillation_fitter.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/jake/Projects/Fitter/StatOnly/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable oscillation_fitter"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/oscillation_fitter.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/oscillation_fitter.dir/build: oscillation_fitter
.PHONY : CMakeFiles/oscillation_fitter.dir/build

CMakeFiles/oscillation_fitter.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/oscillation_fitter.dir/cmake_clean.cmake
.PHONY : CMakeFiles/oscillation_fitter.dir/clean

CMakeFiles/oscillation_fitter.dir/depend:
	cd /home/jake/Projects/Fitter/StatOnly/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jake/Projects/Fitter/StatOnly /home/jake/Projects/Fitter/StatOnly /home/jake/Projects/Fitter/StatOnly/build /home/jake/Projects/Fitter/StatOnly/build /home/jake/Projects/Fitter/StatOnly/build/CMakeFiles/oscillation_fitter.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/oscillation_fitter.dir/depend


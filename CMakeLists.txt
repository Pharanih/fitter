cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

# Define the project
project(oscillation_fitter VERSION 1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Locate ROOT
find_package(ROOT REQUIRED COMPONENTS RIO Hist Tree MathCore Minuit)
if (NOT ROOT_FOUND)
    message(FATAL_ERROR "ROOT not found. Please install ROOT and make sure it is in your PATH.")
endif()

# Locate Prob3plusplus
set(Prob3plusplus_DIR "/home/jake/Projects/Fitter/StatOnly/Prob3plusplus/build/lib/cmake/Prob3plusplus")
find_package(Prob3plusplus REQUIRED)

if (NOT Prob3plusplus_FOUND)
    message(FATAL_ERROR "Prob3plusplus not found. Please ensure Prob3plusplus is built and installed correctly.")
endif()

# Add include directories
include_directories(
    ${ROOT_INCLUDE_DIRS}
    ${Prob3plusplus_INCLUDE_DIRS}
)

# Add source files
set(SOURCES
    oscillation_fitter.cxx
)

# Add headers
set(HEADERS
    oscillation_fitter.h
)

# Create executable
add_executable(oscillation_fitter ${SOURCES} ${HEADERS})

# Link libraries
target_link_libraries(oscillation_fitter
    ${ROOT_LIBRARIES}
    Prob3plusplus
)

# Specify library search path explicitly
target_link_directories(oscillation_fitter PRIVATE
    /home/jake/Projects/Fitter/StatOnly/Prob3plusplus/build
)

# Ensure runtime linker finds the Prob3plusplus shared library
set_target_properties(oscillation_fitter PROPERTIES
    INSTALL_RPATH "/home/jake/Projects/Fitter/StatOnly/Prob3plusplus/build"
    BUILD_RPATH "/home/jake/Projects/Fitter/StatOnly/Prob3plusplus/build"
)

# Installation
install(TARGETS oscillation_fitter
    RUNTIME DESTINATION bin
)



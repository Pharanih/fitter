#!/bin/bash

# Exit on any error
set -e

# Directories
ROOT_DIR=$(pwd)
BUILD_DIR="$ROOT_DIR/build"
PROB3PLUSPLUS_DIR="$ROOT_DIR/Prob3plusplus"
INSTALL_DIR="$PROB3PLUSPLUS_DIR/build"

# Build and install Prob3plusplus
echo "Checking for Prob3plusplus build and installation..."
if [ -f "$INSTALL_DIR/lib/libProb3plusplus.so" ]; then
    echo "Prob3plusplus is already built and installed. Skipping this step."
else
    echo "Building and installing Prob3plusplus..."
    cd "$PROB3PLUSPLUS_DIR"
    mkdir -p build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make -j$(nproc)
    make install
    echo "Prob3plusplus installed successfully in $INSTALL_DIR."
fi

# Build the oscillation fitter
echo "Checking for oscillation fitter build..."
if [ -f "$BUILD_DIR/oscillation_fitter" ]; then
    echo "Oscillation fitter is already built. Skipping this step."
else
    echo "Building oscillation fitter..."
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"
    cmake .. -DCMAKE_BUILD_TYPE=Release -DProb3plusplus_DIR="$INSTALL_DIR/lib/cmake/Prob3plusplus"
    make -j$(nproc)
    echo "Oscillation fitter built successfully."
fi

# Run the fitter
echo "Running oscillation fitter..."
cd "$BUILD_DIR"
./oscillation_fitter






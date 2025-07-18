#!/bin/bash

set -e

# --- Build Directories ---
ROOT_DIR=$(pwd)
BUILD_DIR="$ROOT_DIR/build"
CPU_BUILD_DIR="$BUILD_DIR/cpu"
GPU_BUILD_DIR="$BUILD_DIR/gpu"

# --- Clean and Create Directories ---
echo "--- Cleaning up old builds and creating directories ---"
rm -rf "$BUILD_DIR"
mkdir -p "$CPU_BUILD_DIR"
mkdir -p "$GPU_BUILD_DIR"
echo "Done."
echo ""

# --- Build CPU Version ---
echo "--- Building CPU Version ---"
cd "$CPU_BUILD_DIR"
cmake -DBUILD_GPU=OFF "$ROOT_DIR"
make
echo "CPU build complete."
mv ./cpu_kmer_counter "$ROOT_DIR/"
echo "CPU executable moved to $ROOT_DIR"
echo ""

# --- Build GPU Version ---
echo "--- Building GPU Version ---"
cd "$GPU_BUILD_DIR"
cmake -DBUILD_GPU=ON "$ROOT_DIR"
make
echo "GPU build complete."
mv ./gpu_kmer_counter "$ROOT_DIR/"
echo "GPU executable moved to $ROOT_DIR"
echo ""

echo "----------------------------------------------------"
echo "All builds finished successfully!"
echo "Executables are located in the '$ROOT_DIR' directory."
echo "----------------------------------------------------"

cd "$ROOT_DIR"
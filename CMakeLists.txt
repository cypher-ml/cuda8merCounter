# Set the minimum required version of CMake.
cmake_minimum_required(VERSION 3.18)

# --- Project Definition ---
# Define the project and enable CXX by default. CUDA will be enabled conditionally later.
project(kmer_counter LANGUAGES CXX)

# --- General Project Settings ---
# These settings are common to both CPU and GPU builds.
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# High-performance flags for release builds.
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -mtune=native -flto -ffast-math")

# Add the 'include' directory for header files.
include_directories(include)

# Find and require the pthreads library. This now works because a language is enabled.
find_package(Threads REQUIRED)


# --- Build Target Selection ---
# This creates a user-configurable option.
# Usage: cmake .. -DBUILD_GPU=ON  (to enable)
#        cmake .. -DBUILD_GPU=OFF (to disable, default)
option(BUILD_GPU "Enable to build the GPU-accelerated version" OFF)


# --- Define Source File Sets ---
# Common utilities used by both CPU and GPU versions.
set(COMMON_SOURCES
    src/utils/fastaParser.cpp
    src/utils/out.cpp
)

# Source files specific to the original CPU implementation.
set(CPU_SOURCES
    src/cpu/mainCpu.cpp
    src/cpu/cpuMultithread.cpp
    src/cpu/cpuPipeline.cpp
    src/cpu/kmerCounter.cpp
    ${COMMON_SOURCES}
)

# Source files for the new GPU implementation.
# We need to include cpuMultithread.cpp because the GPU version reuses the producer logic.
set(GPU_SOURCES
    src/gpu/mainGpu.cu
    src/cpu/cpuMultithread.cpp
    ${COMMON_SOURCES}
)


# --- Conditional Build Configuration ---
# This block checks the value of the BUILD_GPU option.
if(BUILD_GPU)
    # --- GPU Build Configuration ---
    message(STATUS "Configuring for GPU build...")

    # Enable the CUDA language for the project.
    enable_language(CUDA)
    
    # Manually specify the C++ standard for the NVCC compiler.
    # This is more robust than relying on CMAKE_CUDA_STANDARD.
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --std=c++20")

    # Specify target GPU architectures for compilation.
    # Adjust these to match your specific GPU for best performance.
    set(CMAKE_CUDA_ARCHITECTURES "75;80;86;90")

    # Define the executable for the GPU version.
    add_executable(gpu_kmer_counter ${GPU_SOURCES})
    
    # Link necessary libraries.
    target_link_libraries(gpu_kmer_counter PRIVATE Threads::Threads)

else()
    # --- CPU Build Configuration ---
    message(STATUS "Configuring for CPU-only build...")

    # Define the executable for the CPU version.
    add_executable(cpu_kmer_counter ${CPU_SOURCES})
    
    # Link necessary libraries.
    target_link_libraries(cpu_kmer_counter PRIVATE Threads::Threads)

endif()

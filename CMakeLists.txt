cmake_minimum_required(VERSION 3.18)
project(kmer_counter LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -mtune=native -flto -ffast-math")

include_directories(include)
find_package(Threads REQUIRED)

option(BUILD_GPU "Enable to build the GPU-accelerated version" OFF)

set(COMMON_SOURCES
    src/utils/fastaParser.cpp
    src/utils/out.cpp
)

set(CPU_SOURCES
    src/cpu/mainCpu.cpp
    src/cpu/cpuMultithread.cpp
    src/cpu/cpuPipeline.cpp
    src/cpu/kmerCounter.cpp
    ${COMMON_SOURCES}
)

set(GPU_SOURCES
    src/gpu/mainGpu.cu
    ${COMMON_SOURCES}
)

if(BUILD_GPU)
    message(STATUS "Configuring for GPU build...")
    enable_language(CUDA)
    
    # CUDA 12.9 has good C++20 support
    set(CMAKE_CUDA_STANDARD 20)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    
    # Updated architectures for newer GPUs (you can adjust based on your target hardware)
    set(CMAKE_CUDA_ARCHITECTURES "70;75;80;86;89;90")
    
    # Optimized CUDA flags for CUDA 12.9
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 --use_fast_math --extended-lambda --expt-relaxed-constexpr")
    
    add_executable(gpu_kmer_counter ${GPU_SOURCES})
    target_link_libraries(gpu_kmer_counter PRIVATE Threads::Threads)
    
    # Set properties for mixed C++/CUDA compilation
    set_target_properties(gpu_kmer_counter PROPERTIES
        CXX_STANDARD 20
        CUDA_STANDARD 20
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_RESOLVE_DEVICE_SYMBOLS ON
    )
else()
    message(STATUS "Configuring for CPU-only build...")
    add_executable(cpu_kmer_counter ${CPU_SOURCES})
    target_link_libraries(cpu_kmer_counter PRIVATE Threads::Threads)
endif()
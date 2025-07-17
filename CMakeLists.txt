cmake_minimum_required(VERSION 3.18)

project(kmer_counter LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
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
    project(kmer_counter_gpu LANGUAGES CXX CUDA)
    set(CMAKE_CUDA_STANDARD 20)
    set(CMAKE_CUDA_ARCHITECTURES "70;75;80;86;90")
    message(STATUS "Configuring for GPU build...")
    add_executable(gpu_kmer_counter ${GPU_SOURCES})
    target_link_libraries(gpu_kmer_counter PRIVATE Threads::Threads)

else()
    project(kmer_counter_cpu LANGUAGES CXX)
    message(STATUS "Configuring for CPU-only build...")
    add_executable(cpu_kmer_counter ${CPU_SOURCES})
    target_link_libraries(cpu_kmer_counter PRIVATE Threads::Threads)

endif()
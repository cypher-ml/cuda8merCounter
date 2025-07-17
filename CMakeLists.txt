cmake_minimum_required(VERSION 3.18)
project(kmer_counter LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) 

include_directories(include)

find_package(Threads REQUIRED)

add_executable(cpu_kmer_counter
    src/cpu/mainCpu.cpp
    src/cpu/kmerCounter.cpp
    src/utils/fastaParser.cpp
    src/cpu/cpuPipeline.cpp 
)

target_link_libraries(cpu_kmer_counter PRIVATE Threads::Threads)
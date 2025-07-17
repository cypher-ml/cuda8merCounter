# Set minimum CMake version and project name
cmake_minimum_required(VERSION 3.18)
#project(kmer_counter LANGUAGES CXX CUDA)
project(kmer_counter LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set CUDA standard
#set(CMAKE_CUDA_STANDARD 17)
#set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# Add the 'include' directory for headers
include_directories(include)

# --- Main Application Target ---
# This will be your final program.

#add_executable(kmer_counter
#    src/main.cpp
#    src/utils/fasta_parser.cpp 
#    src/utils/encoder.cpp
#    # Add other .cpp and .cu files here as you create them
#)

# --- Test Target for the Parser ---
# This creates a separate small program just for testing the parser.
add_executable(test_parser
    src/test_parser.cpp          # A main() dedicated to testing
    src/utils/fasta_parser.cpp   # Link the implementation being tested
)
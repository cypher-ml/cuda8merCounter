# Set minimum CMake version and project name.
cmake_minimum_required(VERSION 3.18)
project(kmer_counter LANGUAGES CXX)

# Set the C++ standard to C++23 and require it.
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # Good practice: disable non-standard compiler extensions.

# Add the 'include' directory so source files can find .hpp headers.
# e.g., #include "utils/fastaParser.hpp"
include_directories(include)

# Find the C++ threading library, which is required for std::thread, etc.
find_package(Threads REQUIRED)

# --- Main CPU Application Target ---
# This creates the final cpu_kmer_counter executable by compiling and
# linking all the necessary source files.
add_executable(cpu_kmer_counter
    src/cpu/cpu.cpp
    src/cpu/kmerCounter.cpp
    src/utils/fastaParser.cpp
)

# Link the required threading library to our executable.
target_link_libraries(cpu_kmer_counter PRIVATE Threads::Threads)

# --- (Optional) Test Target ---
# The 'test_parser' executable from your old file is no longer needed
# by the main application, but you can uncomment it if you want to
# keep it for separate, isolated testing.
#
# add_executable(test_parser
#     src/test_parser.cpp
#     src/utils/fastaParser.cpp
# )
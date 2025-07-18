#pragma once

#include "fastaParser.hpp"
#include <string>
#include <vector>
#include <cstdint>

struct CUstream_st; 

void run_gpu_pipeline(const std::string& filepath);

void setup_gpu_resources(
    unsigned long long*& d_histogram, 
    uint8_t** d_buffers, 
    CUstream_st** streams, 
    int num_streams,
    size_t chunk_size_bytes, 
    size_t histogram_size
);

void process_chunks_on_gpu(
    ParallelFastaReader& reader,
    unsigned long long* d_histogram,
    uint8_t** d_buffers,
    CUstream_st** streams,
    int num_streams
);

void display_top_kmers(const std::vector<unsigned long long>& h_histogram);

void cleanup_gpu_resources(
    unsigned long long* d_histogram, 
    uint8_t** d_buffers, 
    CUstream_st** streams, 
    int num_streams
);
#pragma once

#include "utils/constants.hpp"
#include <cstdint>

__global__ void count_kmers_kernel(const uint8_t* d_encoded_data,
                                  size_t total_bases,
                                  unsigned long long* d_histogram,
                                  size_t global_chunk_start = 0,
                                  bool is_first_chunk = true) {
    using namespace FastaUtils;
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = gridDim.x * blockDim.x;
    
    size_t start_pos = 0;
    if (!is_first_chunk) {
        start_pos = K_MER_SIZE - 1;
    }
    
    size_t end_pos = total_bases;
    if (end_pos < K_MER_SIZE) {
        return;
    }
    end_pos = end_pos - (K_MER_SIZE - 1);
    
    for (size_t i = start_pos + idx; i < end_pos; i += stride) {
        uint64_t kmer_bit_buffer = 0;
        bool valid_kmer = true;
        
        for (int j = 0; j < K_MER_SIZE; ++j) {
            size_t base_pos = i + j;
            
            if (base_pos >= total_bases) {
                valid_kmer = false;
                break;
            }
            
            int byte_idx = base_pos / BASES_PER_BYTE;
            int bit_shift = (BASES_PER_BYTE - 1 - (base_pos % BASES_PER_BYTE)) * BITS_PER_BASE;
            uint8_t base_code = (d_encoded_data[byte_idx] >> bit_shift) & BASE_MASK;
            
            kmer_bit_buffer = (kmer_bit_buffer << BITS_PER_BASE) | base_code;
        }
        
        if (valid_kmer) {
            atomicAdd(&d_histogram[kmer_bit_buffer], 1ULL);
        }
    }
}
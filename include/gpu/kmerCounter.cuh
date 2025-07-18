#pragma once

#include "utils/constants.hpp"
#include <cstdint>

__global__ void count_kmers_kernel(const uint8_t* d_encoded_data, 
                                   size_t total_bases, 
                                   unsigned long long* d_histogram) {
    
    using namespace FastaUtils;
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = gridDim.x * blockDim.x;

    for (size_t i = idx; i < total_bases - (K_MER_SIZE - 1); i += stride) {
        uint64_t kmer_bit_buffer = 0;

        for (int j = 0; j < K_MER_SIZE; ++j) {
            size_t base_pos = i + j;
            int byte_idx = base_pos / BASES_PER_BYTE;
            int bit_shift = (BASES_PER_BYTE - 1 - (base_pos % BASES_PER_BYTE)) * BITS_PER_BASE;
            
            uint8_t base_code = (d_encoded_data[byte_idx] >> bit_shift) & BASE_MASK;
            kmer_bit_buffer = (kmer_bit_buffer << BITS_PER_BASE) | base_code;
        }

        atomicAdd(&d_histogram[kmer_bit_buffer], 1ULL);
    }
}

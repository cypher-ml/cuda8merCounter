#pragma once

#include "utils/constants.hpp"
#include <cstdint>

// d_ is a common prefix for device (GPU) pointers
__global__ void count_kmers_kernel(const uint8_t* d_encoded_data, 
                                   size_t total_bases, 
                                   uint64_t* d_histogram) {
    
    using namespace FastaUtils;
    
    // Calculate the global thread ID
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = gridDim.x * blockDim.x;

    // Each thread processes multiple k-mers
    for (size_t i = idx; i < total_bases - (K_MER_SIZE - 1); i += stride) {
        uint64_t kmer_bit_buffer = 0;

        // This part is tricky on GPU. A simpler, though less performant,
        // approach is to re-read the bits for each k-mer.
        // A more optimized version would use a sliding window in registers.
        for (int j = 0; j < K_MER_SIZE; ++j) {
            size_t base_pos = i + j;
            int byte_idx = base_pos / BASES_PER_BYTE;
            int bit_shift = (BASES_PER_BYTE - 1 - (base_pos % BASES_PER_BYTE)) * BITS_PER_BASE;
            
            uint8_t base_code = (d_encoded_data[byte_idx] >> bit_shift) & BASE_MASK;
            kmer_bit_buffer = (kmer_bit_buffer << BITS_PER_BASE) | base_code;
        }

        // Use atomicAdd to safely increment the counter for this k-mer
        atomicAdd(&d_histogram[kmer_bit_buffer], 1);
    }
}
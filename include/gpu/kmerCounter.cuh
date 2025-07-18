#pragma once

#include "utils/constants.hpp"
#include <cstdint>

/**
 * @brief CUDA kernel to count 8-mers in parallel on the GPU.
 * @details Each thread in the grid processes a subset of the possible k-mer starting
 * positions in the sequence. It constructs the 8-mer, calculates its
 * corresponding index, and uses an atomic operation (atomicAdd) to safely
 * increment the count in the global histogram. The kernel handles chunk
 * boundaries by adjusting its starting position to avoid double-counting
 * k-mers from the overlap region.
 *
 * @param d_encoded_data Pointer to the encoded sequence data in GPU global memory.
 * @param total_bases The total number of valid bases in the sequence data.
 * @param d_histogram Pointer to the global histogram in GPU memory where counts are stored.
 * @param is_first_chunk A boolean flag. If false, the kernel skips the first K-1 bases
 * to avoid recounting k-mers from the chunk overlap.
 */
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
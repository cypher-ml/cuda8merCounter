#include "cpu/kmerCounter.hpp"
#include "utils/constants.hpp"


void count_kmers_in_chunk_with_boundaries(const EncodedChunk& chunk, 
                                         Histogram& histogram,
                                         bool is_first_chunk) {
    if (chunk.data.empty() || chunk.base_count < FastaUtils::K_MER_SIZE) {
        return;
    }

    const auto& data = chunk.data;
    uint64_t kmer_bit_buffer = 0;
    
    size_t start_offset = 0;
    if (!is_first_chunk && chunk.base_count > FastaUtils::K_MER_SIZE - 1) {
        // Skip the first K_MER_SIZE-1 k-mers to avoid double counting
        start_offset = FastaUtils::K_MER_SIZE - 1;
    }

    // Build initial k-mer
    for (int i = 0; i < FastaUtils::K_MER_SIZE - 1; ++i) {
        int byte_idx = i / FastaUtils::BASES_PER_BYTE;
        int bit_shift = (FastaUtils::BASES_PER_BYTE - 1 - (i % FastaUtils::BASES_PER_BYTE)) * FastaUtils::BITS_PER_BASE;
        uint8_t base_code = (data[byte_idx] >> bit_shift) & FastaUtils::BASE_MASK;
        kmer_bit_buffer = (kmer_bit_buffer << 2) | base_code;
    }

    // Process k-mers, skipping overlap region if not first chunk
    for (size_t i = FastaUtils::K_MER_SIZE - 1; i < chunk.base_count; ++i) {
        int byte_idx = i / FastaUtils::BASES_PER_BYTE;
        int bit_shift = (FastaUtils::BASES_PER_BYTE - 1 - (i % FastaUtils::BASES_PER_BYTE)) * FastaUtils::BITS_PER_BASE;
        uint8_t base_code = (data[byte_idx] >> bit_shift) & FastaUtils::BASE_MASK;

        kmer_bit_buffer = ((kmer_bit_buffer << 2) | base_code) & FastaUtils::KMER_MASK;

        // Only count k-mers that start at or after the start_offset
        if (i >= start_offset + FastaUtils::K_MER_SIZE - 1) {
            histogram[kmer_bit_buffer]++;
        }
    }
}
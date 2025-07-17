#include "cpu/kmerCounter.hpp"

constexpr uint64_t KMER_MASK = (1ULL << (K_MER_SIZE * 2)) - 1;

void count_kmers_in_chunk(const EncodedChunk& chunk, Histogram& histogram) {
    if (chunk.data.empty() || chunk.base_count < K_MER_SIZE) {
        return;
    }

    const auto& data = chunk.data;
    uint64_t kmer_bit_buffer = 0;

    for (int i = 0; i < K_MER_SIZE - 1; ++i) {
        int byte_idx = i / FastaUtils::BASES_PER_BYTE;
        int bit_shift = (FastaUtils::BASES_PER_BYTE - 1 - (i % FastaUtils::BASES_PER_BYTE)) * FastaUtils::BITS_PER_BASE;
        uint8_t base_code = (data[byte_idx] >> bit_shift) & FastaUtils::BASE_MASK;
        kmer_bit_buffer = (kmer_bit_buffer << 2) | base_code;
    }

    for (size_t i = K_MER_SIZE - 1; i < chunk.base_count; ++i) {
        int byte_idx = i / FastaUtils::BASES_PER_BYTE;
        int bit_shift = (FastaUtils::BASES_PER_BYTE - 1 - (i % FastaUtils::BASES_PER_BYTE)) * FastaUtils::BITS_PER_BASE;
        uint8_t base_code = (data[byte_idx] >> bit_shift) & FastaUtils::BASE_MASK;
        
        kmer_bit_buffer = ((kmer_bit_buffer << 2) | base_code) & KMER_MASK;

        histogram[kmer_bit_buffer]++;
    }
}
#include "cpu/kmerCounter.hpp"

// We need a mask to extract a single 8-mer (8 bases * 2 bits/base = 16 bits).
constexpr uint64_t KMER_MASK = (1ULL << (K_MER_SIZE * 2)) - 1;

void count_kmers_in_chunk(const EncodedChunk& chunk, Histogram& histogram) {
    if (chunk.data.empty() || chunk.base_count < K_MER_SIZE) {
        return;
    }

    const auto& data = chunk.data;
    uint64_t kmer_bit_buffer = 0;

    // 1. Prime the buffer with the first K-1 bases.
    // This requires careful reading across the first few bytes.
    for (int i = 0; i < K_MER_SIZE - 1; ++i) {
        int byte_idx = i / FastaUtils::BASES_PER_BYTE;
        int bit_shift = (FastaUtils::BASES_PER_BYTE - 1 - (i % FastaUtils::BASES_PER_BYTE)) * FastaUtils::BITS_PER_BASE;
        uint8_t base_code = (data[byte_idx] >> bit_shift) & FastaUtils::BASE_MASK;
        kmer_bit_buffer = (kmer_bit_buffer << 2) | base_code;
    }

    // 2. Slide the window across the rest of the sequence.
    for (size_t i = K_MER_SIZE - 1; i < chunk.base_count; ++i) {
        int byte_idx = i / FastaUtils::BASES_PER_BYTE;
        int bit_shift = (FastaUtils::BASES_PER_BYTE - 1 - (i % FastaUtils::BASES_PER_BYTE)) * FastaUtils::BITS_PER_BASE;
        uint8_t base_code = (data[byte_idx] >> bit_shift) & FastaUtils::BASE_MASK;
        
        // Add the new base and apply the mask to keep it 16 bits.
        kmer_bit_buffer = ((kmer_bit_buffer << 2) | base_code) & KMER_MASK;

        // The k-mer itself is the index into our histogram.
        histogram[kmer_bit_buffer]++;
    }
}
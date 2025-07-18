#pragma once

#include <cstdint>
#include <array>
#include <vector>
#include <string>

/**
 * @namespace FastaUtils
 * @brief Provides constants and utility functions for FASTA sequence processing.
 * @details This namespace centralizes all the core definitions required for
 * k-mer analysis, including k-mer size, encoding/decoding logic,
 * and data structure definitions.
 */
namespace FastaUtils {
    using Histogram = std::vector<uint64_t>;
    constexpr int K_MER_SIZE = 8;
    constexpr uint8_t INVALID_BASE_CODE = 255;
    constexpr int BITS_PER_BASE = 2;
    constexpr int BASES_PER_BYTE = 4;
    constexpr uint8_t BASE_MASK = 0b11;
    constexpr uint64_t KMER_MASK = (1ULL << (K_MER_SIZE * 2)) - 1;

    /**
     * @brief Creates a lookup table for encoding ASCII DNA characters to a 2-bit representation.
     * @details This function is executed at compile time (`consteval`) to generate a
     * mapping from characters ('A', 'C', 'G', 'T', case-insensitive) to their
     * corresponding 2-bit codes. All other characters are mapped to `INVALID_BASE_CODE`.
     * @return A std::array of 256 elements serving as the encoding lookup table.
     */
    consteval std::array<uint8_t, 256> create_encoding_lut() {
        std::array<uint8_t, 256> table{};
        table.fill(INVALID_BASE_CODE);
        table[static_cast<unsigned char>('A')] = table[static_cast<unsigned char>('a')] = 0b00;
        table[static_cast<unsigned char>('C')] = table[static_cast<unsigned char>('c')] = 0b01;
        table[static_cast<unsigned char>('G')] = table[static_cast<unsigned char>('g')] = 0b10;
        table[static_cast<unsigned char>('T')] = table[static_cast<unsigned char>('t')] = 0b11;
        return table;
    }

    /**
     * @brief Decodes a 16-bit k-mer index back into its DNA string representation.
     * @param kmer_index The 16-bit integer representing the k-mer.
     * @return The decoded 8-character DNA string (e.g., "GATTACA").
     */
    std::string decode_kmer(const uint16_t kmer_index);
    
    /**
     * @brief Saves the final k-mer frequency histogram to a tab-separated values (TSV) file.
     * @param histogram The histogram containing the k-mer counts.
     * @param filepath The path to the output TSV file.
     */
    void save_histogram_to_tsv(const Histogram& histogram, const std::string& filepath);
}

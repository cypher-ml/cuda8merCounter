#pragma once

#include <cstdint>
#include <array>
#include <vector>
#include <string>


namespace FastaUtils {
    using Histogram = std::vector<uint64_t>;
    constexpr int K_MER_SIZE = 8;
    constexpr uint8_t INVALID_BASE_CODE = 255;
    constexpr int BITS_PER_BASE = 2;
    constexpr int BASES_PER_BYTE = 4;
    constexpr uint8_t BASE_MASK = 0b11;
    constexpr uint64_t KMER_MASK = (1ULL << (K_MER_SIZE * 2)) - 1;

    consteval std::array<uint8_t, 256> create_encoding_lut() {
        std::array<uint8_t, 256> table{};
        table.fill(INVALID_BASE_CODE);
        table[static_cast<unsigned char>('A')] = table[static_cast<unsigned char>('a')] = 0b00;
        table[static_cast<unsigned char>('C')] = table[static_cast<unsigned char>('c')] = 0b01;
        table[static_cast<unsigned char>('G')] = table[static_cast<unsigned char>('g')] = 0b10;
        table[static_cast<unsigned char>('T')] = table[static_cast<unsigned char>('t')] = 0b11;
        return table;
    }

    std::string decode_kmer(const uint16_t kmer_index);
    void save_histogram_to_tsv(const Histogram& histogram, const std::string& filepath);
}

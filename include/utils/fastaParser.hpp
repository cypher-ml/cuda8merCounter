#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <cstdint>
#include <array>

// DO NOT put "using namespace std;" here

// A struct to hold a chunk of encoded data and its actual base count.
struct EncodedChunk {
    std::vector<uint8_t> data; // <--- Changed
    size_t base_count;
};

// --- K-mer Constants & Utilities ---
constexpr int K_MER_SIZE = 8;

namespace FastaUtils {
    constexpr uint8_t INVALID_BASE_CODE = 255;
    constexpr int BITS_PER_BASE = 2;
    constexpr int BASES_PER_BYTE = 4;
    constexpr uint8_t BASE_MASK = 0b11;

    consteval std::array<uint8_t, 256> create_encoding_lut() { // <--- Changed
        std::array<uint8_t, 256> table{}; // <--- Changed
        table.fill(INVALID_BASE_CODE);
        table[static_cast<unsigned char>('A')] = table[static_cast<unsigned char>('a')] = 0b00;
        table[static_cast<unsigned char>('C')] = table[static_cast<unsigned char>('c')] = 0b01;
        table[static_cast<unsigned char>('G')] = table[static_cast<unsigned char>('g')] = 0b10;
        table[static_cast<unsigned char>('T')] = table[static_cast<unsigned char>('t')] = 0b11;
        return table;
    }

    std::string decode_kmer(uint16_t kmer_index);
}

class FastaStreamReader {
public:
    explicit FastaStreamReader(const std::string& filepath, size_t chunk_size_bases = 4 * 1024 * 1024); // <--- Changed

    EncodedChunk readNextChunk();
    bool isFinished() const;

private:
    static constexpr std::array<uint8_t, 256> ENCODING_LUT = FastaUtils::create_encoding_lut(); // <--- Changed

    std::ifstream m_file; // <--- Changed
    size_t m_chunk_size_bases;
    bool m_is_finished = false;
    std::vector<uint8_t> m_overlap_buffer; // <--- Changed
};
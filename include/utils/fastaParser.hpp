#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <cstdint>
#include <array>


struct EncodedChunk {
    std::vector<uint8_t> data;
    size_t base_count;
};


constexpr int K_MER_SIZE = 8;
using Histogram = std::vector<uint64_t>;


namespace FastaUtils {
    constexpr uint8_t INVALID_BASE_CODE = 255;
    constexpr int BITS_PER_BASE = 2;
    constexpr int BASES_PER_BYTE = 4;
    constexpr uint8_t BASE_MASK = 0b11;

    consteval std::array<uint8_t, 256> create_encoding_lut() {
        std::array<uint8_t, 256> table{}; 
        table.fill(INVALID_BASE_CODE);
        table[static_cast<unsigned char>('A')] = table[static_cast<unsigned char>('a')] = 0b00;
        table[static_cast<unsigned char>('C')] = table[static_cast<unsigned char>('c')] = 0b01;
        table[static_cast<unsigned char>('G')] = table[static_cast<unsigned char>('g')] = 0b10;
        table[static_cast<unsigned char>('T')] = table[static_cast<unsigned char>('t')] = 0b11;
        return table;
    }

    std::string decode_kmer(uint16_t kmer_index);

    void save_histogram_to_tsv(const Histogram& histogram, const std::string& filepath);
}


class FastaStreamReader {
public:
    explicit FastaStreamReader(const std::string& filepath, size_t chunk_size_bases = 16 * 1024 * 1024);

    EncodedChunk read_next_chunk();
    bool is_finished() const;
    size_t get_file_size() const;
    size_t get_bytes_processed() const;

private:
    static constexpr std::array<uint8_t, 256> ENCODING_LUT = FastaUtils::create_encoding_lut(); 

    std::ifstream m_file;           
    size_t m_chunk_size_bases;      
    bool m_is_finished = false;     
    std::vector<uint8_t> m_overlap_buffer; 
    size_t m_file_size = 0;         
    size_t m_current_pos = 0;       
};
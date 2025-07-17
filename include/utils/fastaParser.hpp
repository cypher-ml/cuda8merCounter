#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <cstdint>
#include <array>
#include <atomic>
#include <mutex>

struct EncodedChunk {
    std::vector<uint8_t> data;
    size_t base_count;
    size_t chunk_id;  // Add this to track chunk order for debugging
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

    std::string decode_kmer(const uint16_t kmer_index);
    void save_histogram_to_tsv(const Histogram& histogram, const std::string& filepath);
}


class ParallelFastaReader {
public:
    explicit ParallelFastaReader(const std::string& filepath, size_t chunk_size_bases = 16 * 1024 * 1024);

    // Called by multiple producer threads - thread-safe
    EncodedChunk readNextChunk();
    
    bool isFinished() const;
    size_t getFileSize() const;
    size_t getBytesProcessed() const;

private:
    static constexpr std::array<uint8_t, 256> ENCODING_LUT = FastaUtils::create_encoding_lut();

    std::string m_filepath;
    size_t m_chunk_size_bases;
    size_t m_file_size;
    
    // Atomic counters for thread-safe chunk distribution
    std::atomic<size_t> m_next_chunk_id{0};
    std::atomic<size_t> m_bytes_processed{0};
    std::atomic<size_t> m_total_chunks{0};
    
    // Mutex for any shared state that needs protection
    mutable std::mutex m_state_mutex;
    
    // Helper method to calculate chunk boundaries
    std::pair<size_t, size_t> getChunkBoundaries(size_t chunk_id) const;
};
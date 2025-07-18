#pragma once

#include "constants.hpp"
#include <fstream>
#include <atomic>
#include <mutex>

/**
 * @struct EncodedChunk
 * @brief Represents a chunk of the FASTA file that has been encoded.
 * @details This struct holds the 2-bit encoded sequence data, the number of
 * valid bases in the chunk, and a unique identifier for the chunk.
 */
struct EncodedChunk {
    std::vector<uint8_t> data;
    size_t base_count;
    size_t chunk_id;  
};

/**
 * @class ParallelFastaReader
 * @brief Reads a FASTA file in parallel chunks for high-performance processing.
 * @details This class is designed to be thread-safe, allowing multiple producer
 * threads to read different sections of a large FASTA file simultaneously.
 * It handles file I/O, strips non-standard characters, and encodes
 * bases into a compact 2-bit format. It also manages chunk boundaries
 * to ensure k-mers spanning across chunks are handled correctly by
 * prepending an overlap of K-1 bases.
 */
class ParallelFastaReader {
public:
    explicit ParallelFastaReader(const std::string& filepath, size_t chunk_size_bases = 16 * 1024 * 1024);

    EncodedChunk readNextChunk();
    
    bool isFinished() const;
    size_t getFileSize() const;
    size_t getBytesProcessed() const;

private:
    static constexpr std::array<uint8_t, 256> ENCODING_LUT = FastaUtils::create_encoding_lut();

    std::string m_filepath;
    size_t m_chunk_size_bases;
    size_t m_file_size;
    
    std::atomic<size_t> m_next_chunk_id{0};
    std::atomic<size_t> m_bytes_processed{0};
    std::atomic<size_t> m_total_chunks{0};
    
    mutable std::mutex m_state_mutex;
    
    std::pair<size_t, size_t> getChunkBoundaries(size_t chunk_id) const;
};
#include "utils/fastaParser.hpp"
#include "utils/constants.hpp"
#include <iostream>
#include <stdexcept>

using namespace std;


ParallelFastaReader::ParallelFastaReader(const string& filepath, size_t chunk_size_bases)
    : m_filepath(filepath), m_chunk_size_bases(chunk_size_bases) {
    
    ifstream file(filepath, ios::binary | ios::ate);
    if (!file.is_open()) {
        throw runtime_error("Error: Could not open file " + filepath);
    }
    m_file_size = file.tellg();
    file.close();
    
    size_t chunk_size_bytes = m_chunk_size_bases / FastaUtils::BASES_PER_BYTE;
    m_total_chunks = (m_file_size + chunk_size_bytes - 1) / chunk_size_bytes;
}


bool ParallelFastaReader::isFinished() const {
    return m_next_chunk_id >= m_total_chunks;
}


size_t ParallelFastaReader::getFileSize() const {
    return m_file_size;
}


size_t ParallelFastaReader::getBytesProcessed() const {
    return m_bytes_processed.load();
}


pair<size_t, size_t> ParallelFastaReader::getChunkBoundaries(size_t chunk_id) const {
    size_t chunk_size_bytes = m_chunk_size_bases / FastaUtils::BASES_PER_BYTE;
    size_t start_pos = chunk_id * chunk_size_bytes;
    size_t end_pos = min(start_pos + chunk_size_bytes, m_file_size);
    
    if (chunk_id > 0 && start_pos >= (FastaUtils::K_MER_SIZE - 1)) {
        start_pos -= (FastaUtils::K_MER_SIZE - 1);
    }
    
    return {start_pos, end_pos};
}


EncodedChunk ParallelFastaReader::readNextChunk() {
    size_t chunk_id = m_next_chunk_id.fetch_add(1);
    
    if (chunk_id >= m_total_chunks) {
        return {{}, 0, chunk_id};
    }
    
    size_t chunk_size_bytes = m_chunk_size_bases / FastaUtils::BASES_PER_BYTE;
    size_t start_pos = chunk_id * chunk_size_bytes;
    size_t end_pos = min(start_pos + chunk_size_bytes, m_file_size);
    
    size_t read_start_pos = start_pos;
    bool needs_overlap = false;
    size_t overlap_bases = 0;
    
    if (chunk_id > 0 && start_pos >= (FastaUtils::K_MER_SIZE - 1)) {
        read_start_pos = start_pos - (FastaUtils::K_MER_SIZE - 1);
        needs_overlap = true;
        overlap_bases = FastaUtils::K_MER_SIZE - 1;
    }
    
    ifstream file(m_filepath, ios::binary);
    if (!file.is_open()) {
        throw runtime_error("Error: Could not open file in thread");
    }
    
    file.seekg(read_start_pos);
    
    if (read_start_pos == 0 && file.peek() == '>') {
        string header_line;
        getline(file, header_line);
        read_start_pos = file.tellg();
    }
    
    vector<uint8_t> encoded_data;
    encoded_data.reserve((end_pos - read_start_pos) / FastaUtils::BASES_PER_BYTE + 100);
    
    size_t total_bases_read = 0;
    size_t bases_in_actual_chunk = 0;
    uint8_t packed_byte = 0;
    int bases_in_byte = 0;
    
    const size_t read_buffer_size = 4 * 1024 * 1024;
    vector<char> buffer(read_buffer_size);
    
    size_t current_pos = read_start_pos;
    bool in_overlap_region = needs_overlap;
    
    while (current_pos < end_pos && !file.eof()) {
        size_t to_read = min(read_buffer_size, end_pos - current_pos);
        file.read(buffer.data(), to_read);
        size_t bytes_read = file.gcount();
        
        if (bytes_read == 0) {
            break;
        }
        
        for (size_t i = 0; i < bytes_read; ++i) {
            const char c = buffer[i];
            const uint8_t base_code = ENCODING_LUT[static_cast<unsigned char>(c)];
            
            if (base_code != FastaUtils::INVALID_BASE_CODE) {
                total_bases_read++;
                
                if (in_overlap_region && total_bases_read > overlap_bases) {
                    in_overlap_region = false;
                }
                if (!in_overlap_region) {
                    bases_in_actual_chunk++;
                }
                packed_byte = (packed_byte << FastaUtils::BITS_PER_BASE) | base_code;
                bases_in_byte++;
                
                if (bases_in_byte == FastaUtils::BASES_PER_BYTE) {
                    encoded_data.push_back(packed_byte);
                    packed_byte = 0;
                    bases_in_byte = 0;
                }
            }
        }
        
        current_pos += bytes_read;
    }
    
    if (bases_in_byte > 0) {
        packed_byte <<= (FastaUtils::BASES_PER_BYTE - bases_in_byte) * FastaUtils::BITS_PER_BASE;
        encoded_data.push_back(packed_byte);
    }
    
    m_bytes_processed.fetch_add(end_pos - start_pos);
    
    return {std::move(encoded_data), total_bases_read, chunk_id};
}
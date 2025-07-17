#include "utils/fastaParser.hpp"
#include <iostream>
#include <stdexcept>

// Per your request, this file also uses the std namespace.
using namespace std;

// --- FastaStreamReader Implementation ---

FastaStreamReader::FastaStreamReader(const string& filepath, size_t chunk_size_bases)
    : m_chunk_size_bases(chunk_size_bases) {
    m_file.open(filepath, ios::binary);
    if (!m_file.is_open()) {
        throw runtime_error("Error: Could not open file " + filepath);
    }

    if (m_file.peek() == '>') {
        string header_line;
        getline(m_file, header_line);
    }
}

bool FastaStreamReader::isFinished() const {
    return m_is_finished;
}

EncodedChunk FastaStreamReader::readNextChunk() {
    if (m_is_finished) {
        return {{}, 0};
    }

    vector<uint8_t> encoded_data;
    // Reserve memory for the chunk + overlap to reduce reallocations.
    encoded_data.reserve(m_chunk_size_bases / FastaUtils::BASES_PER_BYTE + K_MER_SIZE);

    // CRITICAL: Prepend the overlap from the *previous* chunk.
    // This ensures k-mers that span chunks are processed correctly.
    encoded_data.insert(encoded_data.end(), m_overlap_buffer.begin(), m_overlap_buffer.end());
    
    size_t bases_in_chunk = m_overlap_buffer.size() * FastaUtils::BASES_PER_BYTE;

    uint8_t packed_byte = 0;
    int bases_in_byte = 0;
    char c;

    while (bases_in_chunk < m_chunk_size_bases && m_file.get(c)) {
        const uint8_t base_code = ENCODING_LUT[static_cast<unsigned char>(c)];

        if (base_code != FastaUtils::INVALID_BASE_CODE) {
            bases_in_chunk++;
            
            packed_byte = (packed_byte << FastaUtils::BITS_PER_BASE) | base_code;
            bases_in_byte++;

            if (bases_in_byte == FastaUtils::BASES_PER_BYTE) {
                encoded_data.push_back(packed_byte);
                packed_byte = 0;
                bases_in_byte = 0;
            }
        }
    }

    // Handle the very last bases in the file if they don't fill a whole byte.
    if (bases_in_byte > 0 && m_file.eof()) {
        packed_byte <<= (FastaUtils::BASES_PER_BYTE - bases_in_byte) * FastaUtils::BITS_PER_BASE;
        encoded_data.push_back(packed_byte);
    }

    // If we read nothing and the file is at its end, we are done.
    if (bases_in_chunk == (m_overlap_buffer.size() * FastaUtils::BASES_PER_BYTE) && m_file.eof()) {
        m_is_finished = true;
        return {{}, 0};
    }
    
    // Update the overlap buffer for the *next* call.
    m_overlap_buffer.clear();
    if (encoded_data.size() >= K_MER_SIZE - 1) {
        m_overlap_buffer.assign(encoded_data.end() - (K_MER_SIZE - 1), encoded_data.end());
    }

    if (m_file.eof()) {
        m_is_finished = true;
    }

    return {std::move(encoded_data), bases_in_chunk};
}


std::string FastaUtils::decode_kmer(uint16_t kmer_index) {
    static const std::array<char, 4> DECODING_LUT = {'A', 'C', 'G', 'T'};
    std::string kmer_str(K_MER_SIZE, ' ');
    
    for (int i = K_MER_SIZE - 1; i >= 0; --i) {
        // Get the 2-bit code for the current base.
        uint8_t base_code = kmer_index & BASE_MASK; // BASE_MASK is 0b11
        kmer_str[i] = DECODING_LUT[base_code];
        
        // Move to the next base.
        kmer_index >>= BITS_PER_BASE; // BITS_PER_BASE is 2
    }
    return kmer_str;
}
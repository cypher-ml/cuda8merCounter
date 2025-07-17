#include "utils/fastaParser.hpp"
#include <iostream>
#include <stdexcept>

using namespace std;


FastaStreamReader::FastaStreamReader(const string& filepath, size_t chunk_size_bases)
    : m_chunk_size_bases(chunk_size_bases) {
    m_file.open(filepath, ios::binary | ios::ate);
    if (!m_file.is_open()) {
        throw runtime_error("Error: Could not open file " + filepath);
    }
    m_file_size = m_file.tellg(); 
    m_file.seekg(0, ios::beg);   

    if (m_file.peek() == '>') {
        string header_line;
        getline(m_file, header_line);
    }
    m_current_pos = m_file.tellg();
}


bool FastaStreamReader::is_finished() const {
    return m_is_finished;
}


size_t FastaStreamReader::get_file_size() const {
    return m_file_size;
}


size_t FastaStreamReader::get_bytes_processed() const {
    return m_current_pos;
}


EncodedChunk FastaStreamReader::read_next_chunk() {
    if (m_is_finished) {
        return {{}, 0};
    }

    vector<uint8_t> encoded_data;
    encoded_data.reserve(m_chunk_size_bases / FastaUtils::BASES_PER_BYTE + K_MER_SIZE);
    encoded_data.insert(encoded_data.end(), m_overlap_buffer.begin(), m_overlap_buffer.end());
    
    size_t bases_in_chunk = m_overlap_buffer.size() * FastaUtils::BASES_PER_BYTE;
    uint8_t packed_byte = 0;
    int bases_in_byte = 0;

    const size_t read_buffer_size = 1024 * 1024;
    vector<char> buffer(read_buffer_size);

    while (bases_in_chunk < m_chunk_size_bases && !m_file.eof()) {
        m_file.read(buffer.data(), buffer.size());
        size_t bytes_read = m_file.gcount();

        if (bytes_read == 0) {
            break;
        }
        for (size_t i = 0; i < bytes_read; ++i) {
            const char c = buffer[i];
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
    }

    m_current_pos = m_file.tellg();
    
    if (bases_in_byte > 0 && m_file.eof()) {
        packed_byte <<= (FastaUtils::BASES_PER_BYTE - bases_in_byte) * FastaUtils::BITS_PER_BASE;
        encoded_data.push_back(packed_byte);
    }

    if (bases_in_chunk == (m_overlap_buffer.size() * FastaUtils::BASES_PER_BYTE) && m_file.eof()) {
        m_is_finished = true;
        return {{}, 0};
    }
    
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
        uint8_t base_code = kmer_index & BASE_MASK; 
        kmer_str[i] = DECODING_LUT[base_code];
        
        kmer_index >>= BITS_PER_BASE; 
    }
    return kmer_str;
}


void FastaUtils::save_histogram_to_tsv(const Histogram& histogram, const std::string& filepath) {
    cout << "\nSaving full histogram to " << filepath << "..." << endl;
    ofstream out_file(filepath);
    if (!out_file) {
        cerr << "Error: Could not open file for writing: " << filepath << endl;
        return;
    }

    out_file << "k-mer\tcount\n";

    for (size_t i = 0; i < histogram.size(); ++i) {
        if (histogram[i] > 0) {
            out_file << FastaUtils::decode_kmer(i) << "\t" << histogram[i] << "\n";
        }
    }
    cout << "Finished saving histogram." << endl;
}
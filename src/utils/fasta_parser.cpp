#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <array>

using namespace std;


consteval array<char, 256> create_lut() {
    array<char, 256> lut{};
    lut[static_cast<unsigned char>('A')] = 'A'; lut[static_cast<unsigned char>('C')] = 'C';
    lut[static_cast<unsigned char>('G')] = 'G'; lut[static_cast<unsigned char>('T')] = 'T';
    lut[static_cast<unsigned char>('a')] = 'A'; lut[static_cast<unsigned char>('c')] = 'C';
    lut[static_cast<unsigned char>('g')] = 'G'; lut[static_cast<unsigned char>('t')] = 'T';
    return lut;
}

string read_and_clean_fasta(const string& filepath) {
    ifstream file(filepath, ios::binary);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filepath << endl;
        return "";
    }

    if (file.peek() == '>') {
        string header_line;
        getline(file, header_line);
    }

    streampos sequence_start_pos = file.tellg();

    file.seekg(0, ios::end);
    streampos file_end_pos = file.tellg();
    
    file.seekg(sequence_start_pos);

    string clean_sequence;
    if (file_end_pos > sequence_start_pos) {
        clean_sequence.reserve(static_cast<size_t>(file_end_pos - sequence_start_pos));
    }
    
    static constexpr array<char, 256> lut = create_lut();

    const size_t buffer_size = 1024 * 1024; 
    vector<char> buffer(buffer_size);

    while (file.read(buffer.data(), buffer.size()) || file.gcount() > 0) {
        const char* buf_ptr = buffer.data();
        const size_t bytes_read = file.gcount();
        
        for (size_t i = 0; i < bytes_read; ++i) {
            if (const char mapped_char = lut[static_cast<unsigned char>(buf_ptr[i])]) {
                clean_sequence.push_back(mapped_char);
            }
        }
    }

    return clean_sequence;
}
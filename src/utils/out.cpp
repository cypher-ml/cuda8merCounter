#include "utils/constants.hpp"
#include <iostream>
#include <fstream>

using namespace std;


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
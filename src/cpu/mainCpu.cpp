#include "cpu/cpuPipeline.hpp"
#include <iostream>

using namespace std;


int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <path_to_fasta_file>" << endl;
        return 1;
    }
    string filepath = argv[1];

    try {
        run_cpu_pipeline(filepath);
    } catch (const exception& e) {
        cerr << "An error occurred during execution: " << e.what() << endl;
        return 1;
    }

    return 0;
}
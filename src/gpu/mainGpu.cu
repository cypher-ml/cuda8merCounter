#include "utils/fastaParser.hpp"
#include "cpu/cpuMultithread.hpp" 
#include "gpu/gpuUtils.cuh"
#include "utils/constants.hpp"

#include <iostream>
#include <chrono>
#include <vector>
#include <thread>
#include <stdexcept>
#include <cuda_runtime.h>

// Bring in your custom using namespace std
using namespace std;

void run_gpu_pipeline(const string& filepath) {
    cout << "Starting GPU-based 8-mer counting on " << filepath << endl;
    auto start_time = chrono::high_resolution_clock::now();

    // === GPU and System Setup ===
    const unsigned int num_producers = thread::hardware_concurrency();
    cout << "Using " << num_producers << " producer threads." << endl;
    const int num_streams = 2;

    // === Device Memory Allocation ===
    unsigned long long* d_histogram = nullptr;
    const size_t histogram_size = (1 << (FastaUtils::K_MER_SIZE * 2));
    uint8_t* d_buffers[num_streams] = {nullptr};
    cudaStream_t streams[num_streams] = {nullptr};
    size_t chunk_size_bytes = (16 * 1024 * 1024);

    setup_gpu_resources(d_histogram, d_buffers, streams, num_streams, chunk_size_bytes, histogram_size);

    // === Producer Thread Setup ===
    ParallelFastaReader reader(filepath);
    vector<thread> producers;
    Threading::active_producers = num_producers;
    Threading::production_finished = false;
    while(!Threading::chunk_queue.empty()) Threading::chunk_queue.pop();

    for (unsigned int i = 0; i < num_producers; ++i) {
        producers.emplace_back(producer_task, ref(reader), i);
    }

    // === Process Chunks ===
    process_chunks_on_gpu(reader, d_histogram, d_buffers, streams, num_streams);

    // Wait for all producers to finish
    for (auto& t : producers) {
        t.join();
    }
    
    // Wait for all GPU work to complete
    cudaDeviceSynchronize();

    auto end_time = chrono::high_resolution_clock::now();
    cout << "Aggregation and computation finished." << endl;

    // === Results Processing ===
    vector<unsigned long long> h_histogram(histogram_size);
    cudaMemcpy(h_histogram.data(), d_histogram, histogram_size * sizeof(unsigned long long), cudaMemcpyDeviceToHost);
    
    auto duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    cout << "Total execution time: " << duration.count() << " seconds." << endl;

    display_top_kmers(h_histogram);
    
    // Save the full histogram to a file
    FastaUtils::save_histogram_to_tsv(
        reinterpret_cast<const FastaUtils::Histogram&>(h_histogram), 
        "kmer_counts_gpu.tsv"
    );

    cleanup_gpu_resources(d_histogram, d_buffers, streams, num_streams);
}


int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <path_to_fasta_file>" << endl;
        return 1;
    }
    string filepath = argv[1];

    try {
        run_gpu_pipeline(filepath);
    } catch (const exception& e) {
        cerr << "An error occurred during execution: " << e.what() << endl;
        return 1;
    }

    return 0;
}
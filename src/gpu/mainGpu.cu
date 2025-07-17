#include "utils/fastaParser.hpp"
#include "cpu/cpuMultithread.hpp" // We'll reuse the producer logic
#include "gpu/kmerCounter.cuh"
#include "utils/constants.hpp"

#include <iostream>
#include <chrono>
#include <vector>
#include <thread>
#include <algorithm>
#include <cuda_runtime.h>

// Bring in your custom using namespace std
using namespace std;

void run_gpu_pipeline(const string& filepath);

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

void run_gpu_pipeline(const string& filepath) {
    cout << "Starting GPU-based 8-mer counting on " << filepath << endl;
    auto start_time = chrono::high_resolution_clock::now();

    // === GPU and System Setup ===
    const unsigned int num_producers = thread::hardware_concurrency();
    cout << "Using " << num_producers << " producer threads." << endl;

    // Define number of concurrent streams for overlapping transfers/kernels
    const int num_streams = 2; 

    // === Device Memory Allocation ===
    unsigned long long* d_histogram; // Changed to unsigned long long
    const size_t histogram_size = (1 << (FastaUtils::K_MER_SIZE * 2));
    cudaMalloc(&d_histogram, histogram_size * sizeof(unsigned long long));
    cudaMemset(d_histogram, 0, histogram_size * sizeof(unsigned long long));

    uint8_t* d_buffers[num_streams];
    cudaStream_t streams[num_streams];
    
    ParallelFastaReader reader(filepath); // The reader for producers
    size_t chunk_size_bytes = (16 * 1024 * 1024); 

    for (int i = 0; i < num_streams; ++i) {
        cudaStreamCreate(&streams[i]);
        // Allocate space for the largest possible chunk
        cudaMalloc(&d_buffers[i], chunk_size_bytes + 1024);
    }

    // === Producer Thread Setup ===
    vector<thread> producers;
    Threading::active_producers = num_producers;
    Threading::production_finished = false;
    while(!Threading::chunk_queue.empty()) Threading::chunk_queue.pop();


    for (unsigned int i = 0; i < num_producers; ++i) {
        producers.emplace_back(producer_task, ref(reader), i);
    }

    int current_stream = 0;
    while (true) {
        unique_lock<mutex> lock(Threading::queue_mutex);
        Threading::cv.wait(lock, [] { return !Threading::chunk_queue.empty() || Threading::production_finished.load(); });

        if (Threading::chunk_queue.empty() && Threading::production_finished.load()) {
            break;
        }

        EncodedChunk chunk = move(Threading::chunk_queue.front());
        Threading::chunk_queue.pop();
        lock.unlock();

        // === Pipeline Stage: Copy and Kernel Launch ===
        // 1. Async Copy Host -> Device
        cudaMemcpyAsync(d_buffers[current_stream], chunk.data.data(), chunk.data.size() * sizeof(uint8_t), cudaMemcpyHostToDevice, streams[current_stream]);

        // 2. Launch Kernel
        int threads_per_block = 256;
        int blocks_per_grid = (chunk.base_count + threads_per_block - 1) / threads_per_block;
        count_kmers_kernel<<<blocks_per_grid, threads_per_block, 0, streams[current_stream]>>>(d_buffers[current_stream], chunk.base_count, d_histogram);
        
        // 3. Switch to the next stream for the next chunk
        current_stream = (current_stream + 1) % num_streams;
    }

    // Wait for all producers to finish
    for (auto& t : producers) {
        t.join();
    }
    
    // Wait for all GPU work to complete
    cudaDeviceSynchronize();

    auto end_time = chrono::high_resolution_clock::now();
    cout << "Aggregation and computation finished." << endl;

    // === Results Processing ===
    vector<unsigned long long> h_histogram(histogram_size); // Changed to unsigned long long
    cudaMemcpy(h_histogram.data(), d_histogram, histogram_size * sizeof(unsigned long long), cudaMemcpyDeviceToHost);
    
    auto duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    cout << "Total execution time: " << duration.count() << " seconds." << endl;

    // Display top 20 (same as your CPU code)
    cout << "\n--- Top 20 Most Frequent 8-mers ---" << endl;
    vector<pair<unsigned long long, uint16_t>> sorted_counts; // Changed to unsigned long long
    for (size_t i = 0; i < h_histogram.size(); ++i) {
        if (h_histogram[i] > 0) {
            sorted_counts.push_back({h_histogram[i], static_cast<uint16_t>(i)});
        }
    }
    sort(sorted_counts.rbegin(), sorted_counts.rend());
    int limit = min(20, (int)sorted_counts.size());
    for (int i = 0; i < limit; ++i) {
        cout << "#" << i + 1 << ":\t"
             << FastaUtils::decode_kmer(sorted_counts[i].second)
             << "\t(Count: " << sorted_counts[i].first << ")" << endl;
    }

    // === Cleanup ===
    cudaFree(d_histogram);
    for(int i = 0; i < num_streams; ++i) {
        cudaStreamDestroy(streams[i]);
        cudaFree(d_buffers[i]);
    }
}

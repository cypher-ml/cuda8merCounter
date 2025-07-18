#include "gpu/gpuUtils.cuh"
#include "gpu/kmerCounter.cuh"
#include "cpu/cpuMultithread.hpp"
#include "utils/constants.hpp"

#include <iostream>
#include <chrono>
#include <vector>
#include <thread>
#include <algorithm>
#include <stdexcept>
#include <cuda_runtime.h>

using namespace std;

void setup_gpu_resources(
    unsigned long long*& d_histogram, 
    uint8_t** d_buffers, 
    cudaStream_t* streams, 
    int num_streams,
    size_t chunk_size_bytes,
    size_t histogram_size
) {
    cout << "Setting up GPU resources..." << endl;
    cudaMalloc(&d_histogram, histogram_size * sizeof(unsigned long long));
    cudaMemset(d_histogram, 0, histogram_size * sizeof(unsigned long long));

    for (int i = 0; i < num_streams; ++i) {
        cudaStreamCreate(&streams[i]);
        cudaMalloc(&d_buffers[i], chunk_size_bytes + 1024); // Allocate space for the largest possible chunk
    }
}


void process_chunks_on_gpu(
    ParallelFastaReader& reader,
    unsigned long long* d_histogram,
    uint8_t** d_buffers,
    cudaStream_t* streams,
    int num_streams
) {
    int current_stream = 0;
    static size_t global_position = 0;  // Track global position across all chunks
    static bool first_chunk = true;
    
    while (true) {
        unique_lock<mutex> lock(Threading::queue_mutex);
        Threading::cv.wait(lock, [] { 
            return !Threading::chunk_queue.empty() || Threading::production_finished.load(); 
        });
        
        if (Threading::chunk_queue.empty() && Threading::production_finished.load()) {
            break;
        }
        
        EncodedChunk chunk = move(Threading::chunk_queue.front());
        Threading::chunk_queue.pop();
        lock.unlock();
        
        // Async Copy Host -> Device
        cudaMemcpyAsync(d_buffers[current_stream], chunk.data.data(), 
                       chunk.data.size() * sizeof(uint8_t), 
                       cudaMemcpyHostToDevice, streams[current_stream]);
        
        // Launch Kernel with boundary information
        int threads_per_block = 256;
        int blocks_per_grid = (chunk.base_count + threads_per_block - 1) / threads_per_block;
        
        count_kmers_kernel<<<blocks_per_grid, threads_per_block, 0, streams[current_stream]>>>(
            d_buffers[current_stream], 
            chunk.base_count, 
            d_histogram,
            global_position,     // Global position for this chunk
            first_chunk         // Is this the first chunk?
        );
        
        // Update global position and first_chunk flag
        global_position += chunk.base_count;
        first_chunk = false;
        
        current_stream = (current_stream + 1) % num_streams;
    }
}


void display_top_kmers(const vector<unsigned long long>& h_histogram) {
    cout << "\n--- Top 20 Most Frequent 8-mers ---" << endl;
    vector<pair<unsigned long long, uint16_t>> sorted_counts;
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
}


void cleanup_gpu_resources(
    unsigned long long* d_histogram, 
    uint8_t** d_buffers, 
    cudaStream_t* streams, 
    int num_streams
) {
    cout << "Cleaning up GPU resources..." << endl;
    cudaFree(d_histogram);
    for(int i = 0; i < num_streams; ++i) {
        cudaStreamDestroy(streams[i]);
        cudaFree(d_buffers[i]);
    }
}
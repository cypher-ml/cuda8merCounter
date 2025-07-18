#pragma once

#include "utils/fastaParser.hpp"
#include <string>
#include <vector>
#include <cstdint>

struct CUstream_st; 

void run_gpu_pipeline(const std::string& filepath);

/**
 * @brief Sets up all necessary GPU resources for the k-mer counting pipeline.
 * @details This includes allocating memory on the device for the main histogram and
 * the data transfer buffers (one for each CUDA stream). It also initializes
 * the CUDA streams used for asynchronous operations.
 * @param d_histogram Reference to the device pointer for the histogram.
 * @param d_buffers Reference to the array of device pointers for chunk data buffers.
 * @param streams Reference to the array of CUDA stream pointers.
 * @param num_streams The number of CUDA streams to create.
 * @param chunk_size_bytes The size of the data buffers to allocate.
 * @param histogram_size The number of elements in the k-mer histogram.
 */
void setup_gpu_resources(
    unsigned long long*& d_histogram, 
    uint8_t** d_buffers, 
    CUstream_st** streams, 
    int num_streams,
    size_t chunk_size_bytes, 
    size_t histogram_size
);

/**
 * @brief Manages the core GPU processing loop.
 * @details This function runs on the host and orchestrates the pipeline:
 * 1. Fetches an encoded chunk from the concurrent queue.
 * 2. Asynchronously copies the chunk to a device buffer using a CUDA stream.
 * 3. Launches the `count_kmers_kernel` in the same stream to process the data.
 * It cycles through the available streams to overlap data transfer with kernel execution.
 * @param reader A reference to the ParallelFastaReader (used for signaling).
 * @param d_histogram Device pointer to the global k-mer histogram.
 * @param d_buffers Array of device pointers to the data buffers.
 * @param streams Array of CUDA streams.
 * @param num_streams The number of CUDA streams being used.
 */
void process_chunks_on_gpu(
    ParallelFastaReader& reader,
    unsigned long long* d_histogram,
    uint8_t** d_buffers,
    CUstream_st** streams,
    int num_streams
);

/**
 * @brief Displays the top N most frequent k-mers from the final histogram.
 * @param h_histogram A reference to the final k-mer histogram on the host.
 */
void display_top_kmers(const std::vector<unsigned long long>& h_histogram);

/**
 * @brief Frees all allocated GPU resources.
 * @details Deallocates the device memory for the histogram and data buffers,
 * and destroys the created CUDA streams.
 * @param d_histogram Device pointer to the histogram.
 * @param d_buffers Array of device pointers to data buffers.
 * @param streams Array of CUDA streams.
 * @param num_streams The number of streams to clean up.
 */
void cleanup_gpu_resources(
    unsigned long long* d_histogram, 
    uint8_t** d_buffers, 
    CUstream_st** streams, 
    int num_streams
);
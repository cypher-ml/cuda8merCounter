#pragma once

#include "utils/fastaParser.hpp" 
#include "cpu/kmerCounter.hpp"
#include <queue>              
#include <mutex>              
#include <condition_variable> 
#include <atomic>             
#include <thread>             
#include <algorithm>          
#include <cstddef>            

/**
 * @namespace Threading
 * @brief Manages shared resources for the multithreaded producer-consumer pipeline.
 * @details This namespace contains the global concurrent queue for encoded chunks,
 * synchronization primitives (mutex, condition variable), and atomic flags
 * to coordinate the producer, consumer, and progress bar threads.
 */
namespace Threading{
    inline std::queue<EncodedChunk> chunk_queue;
    inline std::mutex queue_mutex;
    inline std::condition_variable cv;
    inline std::atomic<bool> production_finished(false);
    inline std::atomic<size_t> active_producers(0);

    inline std::atomic<size_t> g_bytes_processed(0);
    inline std::atomic<size_t> g_file_size(0);

    const unsigned int total_threads = std::thread::hardware_concurrency();
    const unsigned int producer_threads = std::max(2u, total_threads / 2);
    const unsigned int consumer_threads = std::max(1u, total_threads - producer_threads);
}

/**
 * @brief The main function for a producer thread.
 * @details Reads chunks from the ParallelFastaReader and pushes them into the shared queue.
 * @param reader A reference to the ParallelFastaReader instance.
 * @param producer_id A unique identifier for the producer thread (for debugging).
 */
void producer_task(ParallelFastaReader& reader, size_t producer_id);

/**
 * @brief The main function for a consumer thread.
 * @details Fetches encoded chunks from the shared queue and counts the k-mers,
 * storing the results in a local histogram.
 * @param local_histogram A reference to the thread's private histogram for storing k-mer counts.
 */
void consumer_task(Histogram& local_histogram);

/**
 * @brief A task that runs in a separate thread to display a progress bar.
 * @details Periodically checks the number of bytes processed and updates the console.
 * @param reader A reference to the ParallelFastaReader to get progress information.
 */
void progress_bar_task(ParallelFastaReader& reader);
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


void producer_task(ParallelFastaReader& reader, size_t producer_id);
void consumer_task(Histogram& local_histogram);
void progress_bar_task(ParallelFastaReader& reader);
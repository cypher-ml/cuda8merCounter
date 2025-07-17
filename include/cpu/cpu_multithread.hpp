#pragma once

#include "utils/fastaParser.hpp" 
#include "cpu/kmerCounter.hpp" 


namespace Threading{
    queue<EncodedChunk> chunk_queue;
    mutex queue_mutex;
    condition_variable cv;
    atomic<bool> production_finished(false);
    atomic<size_t> active_producers(0);

    atomic<size_t> g_bytes_processed(0);
    atomic<size_t> g_file_size(0);

    const unsigned int total_threads = thread::hardware_concurrency();
    const unsigned int producer_threads = max(2u, total_threads / 2);  // Use half threads as producers
    const unsigned int consumer_threads = total_threads - producer_threads;
}


void producer_task(ParallelFastaReader& reader, size_t producer_id);
void consumer_task(Histogram& local_histogram);
void progress_bar_task(ParallelFastaReader& reader);
#pragma once

#include "utils/fastaParser.hpp"      // For ParallelFastaReader, EncodedChunk, Histogram
#include "cpu/kmerCounter.hpp"        // For count_kmers_in_chunk_with_boundaries

// Always use using namespace std!
using namespace std;


namespace Threading{
    queue<EncodedChunk> chunk_queue;
    mutex queue_mutex;
    condition_variable cv;
    atomic<bool> production_finished(false);
    atomic<size_t> active_producers(0);

    atomic<size_t> g_bytes_processed(0);
    atomic<size_t> g_file_size(0);
}


// Declarations of the multithreaded tasks
void producer_task(ParallelFastaReader& reader, size_t producer_id);
void consumer_task(Histogram& local_histogram);
void progress_bar_task(ParallelFastaReader& reader);
#include "cpu/cpuMultithread.hpp" // This header declares consumer_task
#include "cpu/kmerCounter.hpp"   // This header declares the function it calls

#include <vector>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <atomic>

using namespace std;

// The consumer task is only used by the CPU version. It is now in its own file
// to ensure it is not compiled as part of the GPU build.
void consumer_task(vector<uint64_t>& total_counts) {
    while (true) {
        unique_lock<mutex> lock(Threading::queue_mutex);
        Threading::cv.wait(lock, [] {
            return !Threading::chunk_queue.empty() || Threading::production_finished.load();
        });

        if (Threading::chunk_queue.empty() && Threading::production_finished.load()) {
            break; // All chunks processed
        }

        EncodedChunk chunk = std::move(Threading::chunk_queue.front());
        Threading::chunk_queue.pop();
        lock.unlock();

        // This function is defined in kmerCounter.cpp and is only needed for the CPU consumer.
        count_kmers_in_chunk_with_boundaries(chunk, total_counts, false);
    }
}
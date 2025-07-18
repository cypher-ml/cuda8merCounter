#include "cpu/cpuMultithread.hpp" 
#include "cpu/kmerCounter.hpp"   

#include <vector>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <atomic>

using namespace std;


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

        count_kmers_in_chunk_with_boundaries(chunk, total_counts, chunk.chunk_id == 0);
    }
}
#include "cpu/cpuMultithread.hpp"
#include <iostream>
#include <iomanip>


using namespace std;


void producer_task(ParallelFastaReader& reader, size_t producer_id) {
    while (!reader.isFinished()) {
        EncodedChunk chunk = reader.readNextChunk();
        
        if (!chunk.data.empty()) {
            {
                lock_guard<mutex> lock(Threading::queue_mutex);
                Threading::chunk_queue.push(std::move(chunk));
            }
            Threading::cv.notify_one();
        }
    }
    
    size_t remaining = Threading::active_producers.fetch_sub(1) - 1;
    
    if (remaining == 0) {
        lock_guard<mutex> lock(Threading::queue_mutex);
        Threading::production_finished = true;
        Threading::cv.notify_all();
    }
}

#ifndef __CUDACC__
void consumer_task(Histogram& local_histogram) {
    while (true) {
        EncodedChunk chunk;
        {
            unique_lock<mutex> lock(Threading::queue_mutex);
            Threading::cv.wait(lock, []{ return !Threading::chunk_queue.empty() || Threading::production_finished.load(); });

            if (Threading::chunk_queue.empty() && Threading::production_finished.load()) {
                return;
            }

            chunk = std::move(Threading::chunk_queue.front());
            Threading::chunk_queue.pop();
        }
        
        // Use chunk_id to determine if this is the first chunk
        bool is_first_chunk = (chunk.chunk_id == 0);
        count_kmers_in_chunk_with_boundaries(chunk, local_histogram, is_first_chunk);
    }
}
#endif // __CUDACC__

void progress_bar_task(ParallelFastaReader& reader) {
    Threading::g_file_size = reader.getFileSize();

    while (!Threading::production_finished.load() || Threading::active_producers.load() > 0) {
        size_t file_size = Threading::g_file_size.load();
        size_t bytes_processed = reader.getBytesProcessed();

        if (file_size == 0) {
            this_thread::sleep_for(chrono::milliseconds(50));
            continue;
        }

        float percentage = 100.0f * static_cast<float>(bytes_processed) / static_cast<float>(file_size);
        int bar_width = 50;
        int pos = static_cast<int>(bar_width * percentage / 100.0f);

        cout << "\r[";
        for (int i = 0; i < bar_width; ++i) {
            if (i < pos) cout << "=";
            else if (i == pos) cout << ">";
            else cout << " ";
        }
        cout << "] " << fixed << setprecision(1) << percentage << "% " << flush;

        this_thread::sleep_for(chrono::milliseconds(200));
    }
    cout << "\r[" << string(50, '=') << "] 100.0% " << endl;
}
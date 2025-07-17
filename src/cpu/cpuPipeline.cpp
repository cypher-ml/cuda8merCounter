#include "cpu/cpuPipeline.hpp"
#include "cpu/kmerCounter.hpp"
#include "utils/fastaParser.hpp"

#include <iostream>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <chrono>
#include <vector>
#include <string>
#include <algorithm>
#include <utility>
#include <iomanip>
#include <atomic> 

using namespace std;


namespace {
    queue<EncodedChunk> chunk_queue;
    mutex queue_mutex;
    condition_variable cv;
    atomic<bool> production_finished(false);
    atomic<size_t> active_producers(0);

    atomic<size_t> g_bytes_processed(0);
    atomic<size_t> g_file_size(0);
}


void producer_task(ParallelFastaReader& reader, size_t producer_id) {
    while (!reader.isFinished()) {
        EncodedChunk chunk = reader.readNextChunk();
        
        if (!chunk.data.empty()) {
            {
                lock_guard<mutex> lock(queue_mutex);
                chunk_queue.push(std::move(chunk));
            }
            cv.notify_one();
        }
    }
    
    // Decrement active producers count
    size_t remaining = active_producers.fetch_sub(1) - 1;
    
    // If this was the last producer, signal completion
    if (remaining == 0) {
        lock_guard<mutex> lock(queue_mutex);
        production_finished = true;
        cv.notify_all();
    }
}


void consumer_task(Histogram& local_histogram) {
    while (true) {
        EncodedChunk chunk;
        {
            unique_lock<mutex> lock(queue_mutex);
            cv.wait(lock, []{ return !chunk_queue.empty() || production_finished.load(); });

            if (chunk_queue.empty() && production_finished.load()) {
                return;
            }

            chunk = std::move(chunk_queue.front());
            chunk_queue.pop();
        }
        
        // Use chunk_id to determine if this is the first chunk
        bool is_first_chunk = (chunk.chunk_id == 0);
        count_kmers_in_chunk_with_boundaries(chunk, local_histogram, is_first_chunk);
    }
}


void progress_bar_task(ParallelFastaReader& reader) {
    g_file_size = reader.getFileSize();
    
    while (!production_finished.load() || active_producers.load() > 0) {
        size_t file_size = g_file_size.load();
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



void run_cpu_pipeline(const std::string& filepath) {
    cout << "Starting PARALLEL 8-mer counting on " << filepath << endl;
    auto start_time = chrono::high_resolution_clock::now();

    // Create parallel reader
    ParallelFastaReader reader(filepath);

    const unsigned int total_threads = thread::hardware_concurrency();
    const unsigned int producer_threads = max(2u, total_threads / 2);  // Use half threads as producers
    const unsigned int consumer_threads = total_threads - producer_threads;
    
    cout << "Using " << producer_threads << " producer threads and " 
         << consumer_threads << " consumer threads." << endl;

    vector<thread> producer_threads_vec;
    vector<thread> consumer_threads_vec;
    vector<Histogram> local_histograms(consumer_threads, Histogram(1 << (K_MER_SIZE * 2), 0));

    // Reset global state for this run.
    production_finished = false;
    active_producers = producer_threads;
    g_bytes_processed = 0;
    g_file_size = reader.getFileSize();
    while(!chunk_queue.empty()) chunk_queue.pop();

    // Start progress bar thread
    thread progress_thread(progress_bar_task, ref(reader));

    // Start producer threads
    for (unsigned int i = 0; i < producer_threads; ++i) {
        producer_threads_vec.emplace_back(producer_task, ref(reader), i);
    }

    // Start consumer threads
    for (unsigned int i = 0; i < consumer_threads; ++i) {
        consumer_threads_vec.emplace_back(consumer_task, ref(local_histograms[i]));
    }

    // Wait for all threads to complete
    for (auto& t : producer_threads_vec) {
        t.join();
    }
    
    progress_thread.join();
    
    for (auto& t : consumer_threads_vec) {
        t.join();
    }

    // Aggregate results from all consumer threads into a final histogram.
    cout << "Aggregating results..." << endl;
    Histogram final_histogram(1 << (K_MER_SIZE * 2), 0);
    for (const auto& local_hist : local_histograms) {
        for (size_t i = 0; i < local_hist.size(); ++i) {
            final_histogram[i] += local_hist[i];
        }
    }

    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    cout << "Total execution time: " << duration.count() << " seconds." << endl;

    // Display the top 20 most frequent k-mers.
    cout << "\n--- Top 20 Most Frequent 8-mers ---" << endl;
    vector<pair<uint64_t, uint16_t>> sorted_counts;
    sorted_counts.reserve(final_histogram.size());
    for (size_t i = 0; i < final_histogram.size(); ++i) {
        if (final_histogram[i] > 0) {
            sorted_counts.push_back({final_histogram[i], static_cast<uint16_t>(i)});
        }
    }
    sort(sorted_counts.rbegin(), sorted_counts.rend());

    int limit = min(20, (int)sorted_counts.size());
    for (int i = 0; i < limit; ++i) {
        cout << "#" << i + 1 << ":\t"
             << FastaUtils::decode_kmer(sorted_counts[i].second)
             << "\t(Count: " << sorted_counts[i].first << ")" << endl;
    }

    FastaUtils::save_histogram_to_tsv(final_histogram, "kmer_counts.tsv");
}
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
    bool production_finished = false;

    atomic<size_t> g_bytes_processed(0);
    atomic<size_t> g_file_size(0);
}


void producer_task(FastaStreamReader& reader) {
    g_file_size = reader.getFileSize();

    while (!reader.isFinished()) {
        EncodedChunk chunk = reader.readNextChunk();
        
        g_bytes_processed = reader.getBytesProcessed();

        if (!chunk.data.empty()) {
            {
                lock_guard<mutex> lock(queue_mutex);
                chunk_queue.push(std::move(chunk));
            }
            cv.notify_one();
        }
    }

    {
        lock_guard<mutex> lock(queue_mutex);
        production_finished = true; 
    }
    cv.notify_all();
}


void consumer_task(Histogram& local_histogram) {
    while (true) {
        EncodedChunk chunk;
        {
            unique_lock<mutex> lock(queue_mutex);
            cv.wait(lock, []{ return !chunk_queue.empty() || production_finished; });

            if (chunk_queue.empty() && production_finished) {
                return;
            }

            chunk = std::move(chunk_queue.front());
            chunk_queue.pop();
        }
        count_kmers_in_chunk(chunk, local_histogram);
    }
}


void progress_bar_task(const bool& is_production_finished) {
    while (!is_production_finished) {
        size_t file_size = g_file_size.load();
        size_t bytes_processed = g_bytes_processed.load();

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
    cout << "Starting 8-mer counting on " << filepath << endl;
    auto start_time = chrono::high_resolution_clock::now();

    // The reader is created here and used exclusively by the producer thread.
    FastaStreamReader reader(filepath);

    const unsigned int num_threads = thread::hardware_concurrency();
    cout << "Using " << num_threads << " consumer threads." << endl;

    vector<thread> consumer_threads;
    vector<Histogram> local_histograms(num_threads, Histogram(1 << (K_MER_SIZE * 2), 0));

    // Reset global state for this run.
    production_finished = false;
    g_bytes_processed = 0;
    g_file_size = 0;
    while(!chunk_queue.empty()) chunk_queue.pop();

    // Create and start the producer, consumer, and progress bar threads.
    thread producer(producer_task, ref(reader));
    thread progress_thread(progress_bar_task, ref(production_finished));

    for (unsigned int i = 0; i < num_threads; ++i) {
        consumer_threads.emplace_back(consumer_task, ref(local_histograms[i]));
    }

    // Wait for all threads to complete their work.
    producer.join();
    progress_thread.join();
    for (auto& t : consumer_threads) {
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
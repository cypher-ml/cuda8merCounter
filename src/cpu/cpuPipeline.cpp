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

using namespace std;


namespace {
    queue<EncodedChunk> chunk_queue;
    mutex queue_mutex;
    condition_variable cv;
    bool production_finished = false;
}


void producer_task(const string& filepath) {
    FastaStreamReader reader(filepath);
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


void run_cpu_pipeline(const std::string& filepath) {
    cout << "Starting 8-mer counting on " << filepath << endl;
    auto start_time = chrono::high_resolution_clock::now();

    const unsigned int num_threads = thread::hardware_concurrency();
    cout << "Using " << num_threads << " consumer threads." << endl;

    vector<thread> consumer_threads;
    vector<Histogram> local_histograms(num_threads, Histogram(1 << (K_MER_SIZE * 2), 0));

    thread producer(producer_task, filepath);
    for (unsigned int i = 0; i < num_threads; ++i) {
        consumer_threads.emplace_back(consumer_task, ref(local_histograms[i]));
    }

    producer.join();
    for (auto& t : consumer_threads) {
        t.join();
    }

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

    cout << "\n--- Top 20 Most Frequent 8-mers ---" << endl;
    vector<pair<uint64_t, uint16_t>> sorted_counts;
    sorted_counts.reserve(final_histogram.size());
    for (size_t i = 0; i < final_histogram.size(); ++i) {
        if (final_histogram[i] > 0) {
            sorted_counts.push_back({final_histogram[i], static_cast<uint16_t>(i)});
        }
    }
    sort(sorted_counts.rbegin(), sorted_counts.rend());

    for (size_t i = 0; i < 20 && i < sorted_counts.size(); ++i) {
        cout << "#" << i + 1 << ":\t"
             << FastaUtils::decode_kmer(sorted_counts[i].second)
             << "\t(Count: " << sorted_counts[i].first << ")" << endl;
    }

    FastaUtils::save_histogram_to_tsv(final_histogram, "kmer_counts.tsv");
}
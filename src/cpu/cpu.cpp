#include <iostream>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <chrono>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm> // For std::sort
#include <utility>   // For std::pair

#include "utils/fastaParser.hpp"
#include "cpu/kmerCounter.hpp"

using namespace std;

// ... (producer_task and consumer_task functions remain the same) ...
// --- Shared State for Producer-Consumer ---
queue<EncodedChunk> chunk_queue;
mutex queue_mutex;
condition_variable cv;
bool production_finished = false;

// --- Task Functions ---

void producer_task(const string& filepath) {
    FastaStreamReader reader(filepath);
    while (!reader.isFinished()) {
        EncodedChunk chunk = reader.readNextChunk();
        if (!chunk.data.empty()) {
            {
                lock_guard<mutex> lock(queue_mutex);
                chunk_queue.push(std::move(chunk));
            }
            cv.notify_one(); // Wake up one sleeping consumer.
        }
    }

    { // Signal that production is done.
        lock_guard<mutex> lock(queue_mutex);
        production_finished = true;
    }
    cv.notify_all(); // Wake up all threads so they can exit.
}

void consumer_task(Histogram& local_histogram) {
    while (true) {
        EncodedChunk chunk;
        {
            unique_lock<mutex> lock(queue_mutex);
            // Wait until the queue has something OR production is finished.
            cv.wait(lock, []{ return !chunk_queue.empty() || production_finished; });

            // Exit condition: production is done AND the queue is empty.
            if (chunk_queue.empty() && production_finished) {
                return;
            }

            chunk = std::move(chunk_queue.front());
            chunk_queue.pop();
        } // The lock is released here, allowing other consumers to access the queue.

        // Process the data without holding the lock.
        count_kmers_in_chunk(chunk, local_histogram);
    }
}


// *** NEW HELPER FUNCTION TO SAVE RESULTS ***
void save_histogram_to_tsv(const Histogram& histogram, const string& filepath) {
    cout << "\nSaving full histogram to " << filepath << "..." << endl;
    ofstream out_file(filepath);
    if (!out_file) {
        cerr << "Error: Could not open file for writing: " << filepath << endl;
        return;
    }

    // Write header
    out_file << "k-mer\tcount\n";

    for (size_t i = 0; i < histogram.size(); ++i) {
        if (histogram[i] > 0) { // Optional: only write k-mers that were actually found
            out_file << FastaUtils::decode_kmer(i) << "\t" << histogram[i] << "\n";
        }
    }
    cout << "Finished saving histogram." << endl;
}


int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <path_to_fasta_file>" << endl;
        return 1;
    }
    string filepath = argv[1];

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

    // Aggregate results.
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

    // --- 🔬 REPORTING RESULTS 🔬 ---

    // 1. Find and report the top 20 k-mers.
    cout << "\n--- Top 20 Most Frequent 8-mers ---" << endl;
    vector<pair<uint64_t, uint16_t>> sorted_counts;
    sorted_counts.reserve(final_histogram.size());
    for (size_t i = 0; i < final_histogram.size(); ++i) {
        if (final_histogram[i] > 0) {
            sorted_counts.push_back({final_histogram[i], static_cast<uint16_t>(i)});
        }
    }
    // Sort descending by count.
    sort(sorted_counts.rbegin(), sorted_counts.rend());

    for (size_t i = 0; i < 20 && i < sorted_counts.size(); ++i) {
        cout << "#" << i + 1 << ":\t" 
             << FastaUtils::decode_kmer(sorted_counts[i].second) 
             << "\t(Count: " << sorted_counts[i].first << ")" << endl;
    }

    // 2. Save the full histogram to a TSV file.
    save_histogram_to_tsv(final_histogram, "kmer_counts.tsv");

    return 0;
}
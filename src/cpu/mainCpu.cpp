#include "cpu/kmerCounter.hpp"
#include "utils/fastaParser.hpp"
#include "utils/constants.hpp"
#include "cpu/cpuMultithread.hpp"

#include <iostream>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <utility>
#include <iomanip>

using namespace std;


void run_cpu_pipeline(const std::string& filepath) {
    cout << "Starting PARALLEL 8-mer counting on " << filepath << endl;
    auto start_time = chrono::high_resolution_clock::now();

    // Create parallel reader
    ParallelFastaReader reader(filepath);
    
    cout << "Using " << Threading::producer_threads << " producer threads and " 
         << Threading::consumer_threads << " consumer threads." << endl;

    vector<thread> producer_threads_vec;
    vector<thread> consumer_threads_vec;
    vector<Histogram> local_histograms(Threading::consumer_threads, Histogram(1 << (FastaUtils::K_MER_SIZE * 2), 0));

    // Reset global state for this run.
    Threading::production_finished = false;
    Threading::active_producers = Threading::producer_threads;
    Threading::g_bytes_processed = 0;
    Threading::g_file_size = reader.getFileSize();
    while(!Threading::chunk_queue.empty()) Threading::chunk_queue.pop();

    // Start progress bar thread
    thread progress_thread(progress_bar_task, ref(reader));

    // Start producer threads
    for (unsigned int i = 0; i < Threading::producer_threads; ++i) {
        producer_threads_vec.emplace_back(producer_task, ref(reader), i);
    }

    // Start consumer threads
    for (unsigned int i = 0; i < Threading::consumer_threads; ++i) {
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
    Histogram final_histogram(1 << (FastaUtils::K_MER_SIZE * 2), 0);
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


int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <path_to_fasta_file>" << endl;
        return 1;
    }
    string filepath = argv[1];

    try {
        run_cpu_pipeline(filepath);
    } catch (const exception& e) {
        cerr << "An error occurred during execution: " << e.what() << endl;
        return 1;
    }

    return 0;
}
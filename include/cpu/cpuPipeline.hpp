#pragma once

#include <string>

/**
 * @brief Runs the complete multi-threaded CPU k-mer counting pipeline.
 * * This function orchestrates the producer-consumer model:
 * 1. Starts a producer thread to read and encode the FASTA file in chunks.
 * 2. Starts multiple consumer threads to count k-mers in parallel.
 * 3. Waits for all threads to complete.
 * 4. Aggregates the results from all threads.
 * 5. Reports the top 20 k-mers and saves the full histogram to a file.
 * * @param filepath The path to the input FASTA file.
 */
void run_cpu_pipeline(const std::string& filepath);
#pragma once

#include "utils/fastaParser.hpp"
#include <vector>
#include <cstdint>

// A simple alias for our histogram data structure.
using Histogram = std::vector<uint64_t>;

/**
 * @brief Processes a single encoded chunk and updates a local histogram.
 * @param chunk The EncodedChunk containing packed DNA data.
 * @param histogram The local histogram to update with k-mer counts.
 */
void count_kmers_in_chunk(const EncodedChunk& chunk, Histogram& histogram);
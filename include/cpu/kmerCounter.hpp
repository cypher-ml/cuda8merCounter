#pragma once

#include "utils/fastaParser.hpp"
#include <vector>
#include <cstdint>

using Histogram = std::vector<uint64_t>;

/**
 * @brief Counts all 8-mers within a single encoded chunk.
 * @details This function iterates through the encoded data, constructs 8-mers,
 * and increments their corresponding counts in the histogram. It correctly
 * handles chunk boundaries by skipping the overlapping bases in all but the first chunk
 * to prevent double-counting.
 * @param chunk The EncodedChunk to process.
 * @param histogram The histogram (by reference) where counts will be accumulated.
 * @param is_first_chunk A flag indicating if this is the very first chunk from the file,
 * which does not have a preceding overlap region to skip.
 */
void count_kmers_in_chunk_with_boundaries(const EncodedChunk& chunk, Histogram& histogram, bool is_first_chunk);
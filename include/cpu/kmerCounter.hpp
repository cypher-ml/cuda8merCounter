#pragma once

#include "utils/fastaParser.hpp"
#include <vector>
#include <cstdint>

using Histogram = std::vector<uint64_t>;

void count_kmers_in_chunk_with_boundaries(const EncodedChunk& chunk, Histogram& histogram, bool is_first_chunk);
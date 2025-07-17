#pragma once

#include "utils/fastaParser.hpp"
#include <vector>
#include <cstdint>

using Histogram = std::vector<uint64_t>;

void count_kmers_in_chunk(const EncodedChunk& chunk, Histogram& histogram);
#pragma once
#include "sparse_graph.hpp"

uint32_t reducing_peeling_lower_bound(sparse_graph &g);

bitvector local_search_upper_bound(const sparse_graph &g, size_t iterations = 1);

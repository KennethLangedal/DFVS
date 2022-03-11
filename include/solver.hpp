#pragma once
#include "reductions.hpp"

bitvector solve(graph &g, const bitvector &nodes, graph_search &gs, size_t d);

struct branch_info {
    size_t best_inc = std::numeric_limits<size_t>::max(), best_exc = std::numeric_limits<size_t>::max();
    size_t worst_inc = 0, worst_exc = 0, avg_inc = 0, avg_exc = 0;
};

branch_info analysis(graph &g, graph_search &gs);

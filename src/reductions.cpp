#include "reductions.hpp"
#include <algorithm>
#include <cassert>

graph_search::graph_search(size_t N, bool queue_all)
    : search(num_reductions), visited(num_reductions, bitvector(N)), tmp{} {
    if (!queue_all) {
        for (size_t r = 0; r < num_reductions; ++r) {
            visited[r].fill();
        }
        return;
    }
    for (size_t r = 0; r < num_reductions; ++r) {
        for (size_t u = 0; u < N; ++u) {
            search[r].push_back(u);
        }
    }
}

void push_search(graph_search &gs, size_t u) {
    for (size_t r = 0; r < num_reductions; ++r) {
        if (gs.visited[r].get(u))
            gs.search[r].push_back(u);
        gs.visited[r].reset(u);
    }
}

size_t pop_search(graph_search &gs, size_t rule) {
    assert(rule < num_reductions && !gs.search[rule].empty());
    size_t u = gs.search[rule].back();
    gs.search[rule].pop_back();
    gs.visited[rule].set(u);
    return u;
}
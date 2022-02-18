#pragma once
#include "graph.hpp"

const size_t num_reductions = 6;

enum class reductions { zero_degree, self_edge, one_degree, isolated_vertex, twin_vertices, dominating_vertex };

struct graph_search {
    std::vector<std::vector<size_t>> search;
    std::vector<bitvector<N>> visited;

    graph_search();
};

void add_to_fvs(graph& g, std::vector<size_t>& fvs, graph_search& gs, size_t u);

void deactivate_vertex(graph& g, graph_search& gs, size_t u);

void reduce_graph(graph& g, std::vector<size_t>& fvs, graph_search& gs);
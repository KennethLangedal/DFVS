#pragma once
#include "graph.hpp"

const size_t num_reductions = 10;

enum class reductions { zero_degree,
                        self_edge,
                        one_degree,
                        redundant_edges,
                        isolated_vertex,
                        twin_vertices,
                        dominating_vertex,
                        neighborhood_fold,
                        two_one_fold,
                        specific_pattern };

struct graph_search {
    std::vector<std::vector<size_t>> search;
    std::vector<bitvector<N>> visited;

    bitvector<N> DFS_visited;
    std::vector<size_t> DFS_stack, SCC_id, L;

    graph_search();
};

void add_to_fvs(graph &g, bitvector<N> &fvs, graph_search &gs, size_t u);

void exclude_from_fvs(graph &g, bitvector<N> &fvs, graph_search &gs, size_t u);

void deactivate_vertex(graph &g, graph_search &gs, size_t u);

void reduce_graph(graph &g, bitvector<N> &fvs, const bitvector<N> &nodes, graph_search &gs, std::vector<bitvector<N>> &SCC);
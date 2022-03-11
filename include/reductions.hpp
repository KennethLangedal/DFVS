#pragma once
#include "graph.hpp"

const size_t num_reductions = 9;

enum class reductions { zero_degree,
                        self_edge,
                        one_degree,
                        redundant_edges,
                        isolated_vertex,
                        dominating_vertex,
                        twin_vertices,
                        clique_and_one,
                        square_fold,
                        specific_pattern,
                        redundant_edge_meta };

struct graph_search {
    std::vector<std::vector<size_t>> search;
    std::vector<bitvector> visited;

    bitvector DFS_visited;
    std::vector<size_t> DFS_stack, SCC_id, L;

    bitvector tmp, tmp1, tmp2;

    graph_search(size_t N);
};

void add_to_fvs(graph &g, bitvector &fvs, graph_search &gs, size_t u);

void exclude_from_fvs(graph &g, bitvector &fvs, graph_search &gs, size_t u);

void reduce_graph(graph &g, bitvector &fvs, const bitvector &nodes, graph_search &gs, std::vector<bitvector> &SCC);
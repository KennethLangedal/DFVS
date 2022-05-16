#pragma once
#include "reduction_engine.hpp"

const size_t num_reductions = 11;

enum class reductions { zero_degree,
                        self_edge,
                        one_degree,
                        redundant_edges,
                        isolated_vertex,
                        dominating_vertex,
                        cycle_dominating_vertex,
                        twin_vertices,
                        clique_and_one_fold,
                        square_fold,
                        specific_pattern };

struct graph_search {
    std::vector<std::vector<size_t>> search;
    std::vector<bitvector> visited, SCC;

    bitvector DFS_visited;
    std::vector<size_t> DFS_stack, SCC_id, L;

    std::vector<uint32_t> tmp;

    graph_search(size_t N, bool queue_all = true);
};

void push_search(graph_search &gs, size_t u);

size_t pop_search(graph_search &gs, size_t rule);

void queue_neighbourhood(sparse_graph &g, graph_search &gs, size_t u);

void add_to_fvs(sparse_graph &g, bitvector &fvs, graph_search &gs, reduction_engine &re, size_t u);

void exclude_from_fvs(sparse_graph &g, bitvector &fvs, graph_search &gs, reduction_engine &re, size_t u);

void reduce_graph(sparse_graph &g, bitvector &fvs, graph_search &gs, reduction_engine &re, bool SCC = false);

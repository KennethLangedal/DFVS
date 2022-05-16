#pragma once
#include "sparse_graph.hpp"

enum class reduction {
    add_edge,
    remove_edge,
    remove_vertex,
    fold_clique_and_one,
    fold_square
};

struct action {
    reduction t;
    uint32_t u, v1, v2, w1, w2;
};

class reduction_engine {
private:
    std::vector<action> _log;

    void _undo_add_edge(sparse_graph &g, size_t u, size_t v);
    void _undo_remove_edge(sparse_graph &g, size_t u, size_t v);
    void _undo_remove_vertex(sparse_graph &g, size_t u);

public:
    void add_edge(sparse_graph &g, uint32_t u, uint32_t v);
    void remove_edge(sparse_graph &g, uint32_t u, uint32_t v);
    void remove_include_vertex(sparse_graph &g, uint32_t u);
    void remove_exclude_vertex(sparse_graph &g, uint32_t u);
    void fold_clique_and_one(sparse_graph &g, uint32_t u, uint32_t v, bitvector &fvs);
    void fold_square(sparse_graph &g, uint32_t u, uint32_t v1, uint32_t v2, uint32_t w1, uint32_t w2, bitvector &fvs);

    size_t get_timestamp() const;

    void unfold_graph(sparse_graph &g, size_t time, bitvector &fvs);
};

#pragma once
#include <iostream>
#include <tuple>
#include <vector>

#include "bitvector.hpp"

enum class action {
    add_edge,
    remove_edge,
    remove_vertex,
    fold_clique_and_one,
    fold_square
};

class graph {
private:
    size_t _N;

    std::vector<bitvector> _out_edges, _in_edges;
    bitvector _active;

    std::vector<std::tuple<action, size_t, size_t, size_t, size_t, size_t>> _log;

    void _undo_add_edge(size_t u, size_t v);
    void _undo_remove_edge(size_t u, size_t v);
    void _undo_remove_vertex(size_t u);

public:
    void add_edge(size_t u, size_t v);
    void remove_edge(size_t u, size_t v);
    void remove_vertex(size_t u);
    void fold_clique_and_one(size_t u, size_t v);
    void fold_square(size_t u, size_t v1, size_t v2, size_t w1, size_t w2);

    size_t out_degree(size_t u) const;
    size_t in_degree(size_t u) const;
    bool self_loop(size_t u) const;

    const bitvector &active_vertices() const;
    const bitvector &out(size_t u) const;
    const bitvector &in(size_t u) const;

    size_t get_timestamp() const;
    size_t size() const;

    void unfold_graph(size_t time, bitvector &fvs);

    void print_edgelist(std::ostream &os, const bitvector vertices) const;
    friend std::istream &operator>>(std::istream &is, graph &g);
};

std::istream &operator>>(std::istream &is, graph &g);
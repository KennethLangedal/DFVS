#pragma once
#include <iostream>
#include <tuple>
#include <vector>

#include "bitvector.hpp"

constexpr size_t N = 128;

enum class action {
    add_edge,
    remove_edge,
    remove_vertex,
    fold_neighborhood,
    fold_two_one
};

class graph {
private:
    std::vector<bitvector<N>> out_edges, in_edges;
    bitvector<N> active;

    std::vector<std::tuple<action, size_t, size_t, size_t, size_t>> log;

    void undo_add_edge(size_t u, size_t v);
    void undo_remove_edge(size_t u, size_t v);
    void undo_remove_vertex(size_t u);
    void undo_fold_neighborhood(size_t u, size_t v1, size_t v2);
    void undo_fold_two_one(size_t u, size_t v1, size_t v2, size_t w);

public:
    void add_edge(size_t u, size_t v);
    void remove_edge(size_t u, size_t v);
    void remove_vertex(size_t u);
    void fold_neighborhood(size_t u);
    void fold_two_one(size_t u, size_t v1, size_t v2, size_t w);

    size_t degree(size_t u) const;
    size_t in_degree(size_t u) const;
    bool self_loop(size_t u) const;

    const bitvector<N> &active_vertices() const;
    const bitvector<N> &out(size_t u) const;
    const bitvector<N> &in(size_t u) const;

    size_t get_timestamp() const;
    void unfold_graph(size_t time, bitvector<N> &fvs);

    void print_edgelist(std::ostream &os) const;
    friend std::istream &operator>>(std::istream &is, graph &g);
};

std::istream &operator>>(std::istream &is, graph &g);
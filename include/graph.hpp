#pragma once
#include <iostream>
#include <tuple>
#include <vector>

#include "bitvector.hpp"

constexpr size_t N = 512;

enum class action {
    deactivate_edge,
    deactivate_vertex,
    deactivate_single_out,
    deactivate_single_in,
    fold_neighborhood
};

class graph {
private:
    std::vector<bitvector<N>> out_edges, in_edges;
    bitvector<N> active;

    std::vector<std::tuple<action, size_t, size_t, size_t>> log;

    void undo_deactive_edge(size_t u, size_t v);
    void undo_deactive_vertex(size_t u);
    void undo_deactive_single_out(size_t u, size_t v);
    void undo_deactive_single_in(size_t u, size_t v);
    void undo_fold_neighborhood(size_t u, size_t v1, size_t v2);

public:
    void deactive_edge(size_t u, size_t v);
    void deactive_vertex(size_t u);
    void deactive_single_out(size_t u);
    void deactive_single_in(size_t u);
    void fold_neighborhood(size_t u);

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
#pragma once
#include "bitvector.hpp"
#include <iostream>

class sparse_graph {
private:
    size_t _N;
    bitvector _active;

    std::vector<std::vector<uint32_t>> _out_edges, _out_edges_non_pi, _in_edges, _in_edges_non_pi, _pi_edges;

    std::vector<uint32_t> _original_labels;

public:
    sparse_graph() = default;

    // Relabel constructor
    template <class graph>
    sparse_graph(const graph &g);

    size_t size() const;
    bool is_active(uint32_t u) const;
    bool has_edge(uint32_t u, uint32_t v) const;
    uint32_t original_label(uint32_t u) const;

    const bitvector &active_vertices() const;

    // Warning, only safe from reduction_graph
    void reactivate_vertex(uint32_t u);

    void remove_vertex(uint32_t u);
    void remove_edge(uint32_t u, uint32_t v);
    void add_edge(uint32_t u, uint32_t v);

    size_t out_degree(uint32_t u) const;
    size_t out_degree_non_pi(uint32_t u) const;
    size_t in_degree(uint32_t u) const;
    size_t in_degree_non_pi(uint32_t u) const;
    size_t pi_degree(uint32_t u) const;

    const std::vector<uint32_t> &out(uint32_t u) const;
    const std::vector<uint32_t> &out_non_pi(uint32_t u) const;
    const std::vector<uint32_t> &in(uint32_t u) const;
    const std::vector<uint32_t> &in_non_pi(uint32_t u) const;
    const std::vector<uint32_t> &pi(uint32_t u) const;

    void parse_graph(std::istream &is, size_t N);
};

bool includes(const std::vector<uint32_t> &a, const std::vector<uint32_t> &b);

uint32_t intersection_size(const std::vector<uint32_t> &a, const std::vector<uint32_t> &b);

template <class graph>
sparse_graph::sparse_graph(const graph &g)
    : _N(g.active_vertices().popcount()), _active(_N),
      _out_edges(_N), _out_edges_non_pi(_N), _in_edges(_N), _in_edges_non_pi(_N), _pi_edges(_N),
      _original_labels(g.size()) {

    _active.fill();

    std::vector<uint32_t> new_labels(g.size());
    uint32_t i = 0;
    for (auto u : g.active_vertices()) {
        new_labels[u] = i;
        _original_labels[i++] = u;
    }

    for (auto u : g.active_vertices()) {
        for (auto v : g.out(u))
            _out_edges[new_labels[u]].push_back(new_labels[v]);
        for (auto v : g.out_non_pi(u))
            _out_edges_non_pi[new_labels[u]].push_back(new_labels[v]);
        for (auto v : g.in(u))
            _in_edges[new_labels[u]].push_back(new_labels[v]);
        for (auto v : g.in_non_pi(u))
            _in_edges_non_pi[new_labels[u]].push_back(new_labels[v]);
        for (auto v : g.pi(u))
            _pi_edges[new_labels[u]].push_back(new_labels[v]);
    }
}
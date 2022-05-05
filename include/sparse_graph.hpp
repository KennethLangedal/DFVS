#pragma once
#include "bitvector.hpp"
#include <iostream>
#include <vector>

// Todo, relable constructor

class sparse_graph {
private:
    size_t _N;
    bitvector _active;

    std::vector<std::vector<uint32_t>> _out_edges, _out_edges_non_pi, _in_edges, _in_edges_non_pi, _pi_edges;

public:
    size_t size() const;
    bool is_active(uint32_t u) const;
    bool has_edge(uint32_t u, uint32_t v) const;

    const bitvector &active_vertices() const;

    void remove_include_vertex(uint32_t u);
    void remove_exclude_vertex(uint32_t u);
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

    friend std::istream &operator>>(std::istream &is, sparse_graph &g);
};

std::istream &operator>>(std::istream &is, sparse_graph &g);

bool includes(const std::vector<uint32_t> &a, const std::vector<uint32_t> &b);

uint32_t intersection_size(const std::vector<uint32_t> &a, const std::vector<uint32_t> &b);
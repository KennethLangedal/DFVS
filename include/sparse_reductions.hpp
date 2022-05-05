#pragma once
#include "bitvector.hpp"
#include <cassert>

const size_t num_reductions = 6;

enum class reductions { zero_degree,
                        self_edge,
                        one_degree,
                        redundant_edges,
                        isolated_vertex,
                        dominating_vertex };

struct graph_search {
    std::vector<std::vector<size_t>> search;
    std::vector<bitvector> visited;

    std::vector<uint32_t> tmp;

    graph_search(size_t N, bool queue_all = true);
};

void push_search(graph_search &gs, size_t u);

size_t pop_search(graph_search &gs, size_t rule);

template <class graph>
void queue_neighbourhood(graph &g, graph_search &gs, size_t u) {
    for (auto v : g.out(u))
        push_search(gs, v);
    for (auto v : g.in(u))
        push_search(gs, v);
}

template <class graph>
void add_to_fvs(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    fvs.set(u);
    queue_neighbourhood(g, gs, u);
    g.remove_include_vertex(u);
}

template <class graph>
void exclude_from_fvs(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    queue_neighbourhood(g, gs, u);
    g.remove_exclude_vertex(u);
}

template <class graph>
bool zero_degree_reduction(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    assert(g.is_active(u));
    if (g.in_degree(u) == 0 || g.out_degree(u) == 0) {
        exclude_from_fvs(g, fvs, gs, u);
        return true;
    }
    return false;
}

template <class graph>
bool self_edge_reduction(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    assert(g.is_active(u));
    if (g.has_edge(u, u)) {
        add_to_fvs(g, fvs, gs, u);
        return true;
    }
    return false;
}

template <class graph>
bool one_degree_reduction(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    assert(g.is_active(u));
    if (g.in_degree(u) == 1 || g.out_degree(u) == 1) {
        exclude_from_fvs(g, fvs, gs, u);
        return true;
    }
    return false;
}

template <class graph>
bool redundant_edges(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    assert(g.is_active(u));

    // Zero in or out (excluding 2-cycles)

    if (g.out_degree_non_pi(u) == 0 && g.in_degree_non_pi(u) > 0) {
        gs.tmp.clear();
        for (auto v : g.in_non_pi(u))
            gs.tmp.push_back(v);

        for (auto v : gs.tmp) {
            g.remove_edge(v, u);
            push_search(gs, v);
        }
        push_search(gs, u);
        return true;
    }
    if (g.out_degree_non_pi(u) > 0 && g.in_degree_non_pi(u) == 0) {
        gs.tmp.clear();
        for (auto v : g.out_non_pi(u))
            gs.tmp.push_back(v);

        for (auto v : gs.tmp) {
            g.remove_edge(u, v);
            push_search(gs, v);
        }
        push_search(gs, u);
        return true;
    }

    // DOME rule

    bool res = false;
    gs.tmp.clear();
    for (auto v : g.out_non_pi(u)) {
        if (includes(g.in(v), g.in_non_pi(u)) || includes(g.out(u), g.out_non_pi(v))) {
            gs.tmp.push_back(v);
            push_search(gs, v);
            res = true;
        }
    }

    if (res) {
        for (auto v : gs.tmp) {
            g.remove_edge(u, v);
        }
        push_search(gs, u);
    }

    return res;
}

template <class graph>
bool isolated_vertex_reduction(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    assert(g.is_active(u));
    if (g.pi_degree(u) == g.out_degree(u)) {
        for (auto v : g.pi(u)) {
            if (intersection_size(g.out(u), g.out(v)) != g.out_degree(u) - 1) // Not part of clique
                return false;
        }
        exclude_from_fvs(g, fvs, gs, u);
        return true;
    }
    return false;
}

template <class graph>
bool dominating_vertex_reduction(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    assert(g.is_active(u));

    for (auto v : g.pi(u)) {
        if (intersection_size(g.pi(u), g.out(v)) == g.out_degree(v) - 1 || intersection_size(g.pi(u), g.in(v)) == g.in_degree(v) - 1) {
            add_to_fvs(g, fvs, gs, u);
            return true;
        }
        if (intersection_size(g.pi(v), g.out(u)) == g.out_degree(u) - 1 || intersection_size(g.pi(v), g.in(u)) == g.in_degree(u) - 1) {
            add_to_fvs(g, fvs, gs, v);
            return true;
        }
    }

    return false;
}

template <class graph>
void reduce_graph(graph &g, bitvector &fvs, graph_search &gs) {
    size_t rule = 0;
    while (rule < num_reductions) {
        if (gs.search[rule].empty()) {
            rule++;
        } else {
            size_t u = pop_search(gs, rule);
            if (!g.is_active(u))
                continue;
            bool found = false;
            switch ((reductions)rule) {
            case reductions::zero_degree:
                found = zero_degree_reduction(g, fvs, gs, u);
                break;
            case reductions::self_edge:
                found = self_edge_reduction(g, fvs, gs, u);
                break;
            case reductions::one_degree:
                found = one_degree_reduction(g, fvs, gs, u);
                break;
            case reductions::redundant_edges:
                found = redundant_edges(g, fvs, gs, u);
                break;
            case reductions::isolated_vertex:
                found = isolated_vertex_reduction(g, fvs, gs, u);
                break;
            case reductions::dominating_vertex:
                found = dominating_vertex_reduction(g, fvs, gs, u);
                break;
            default:
                break;
            }
            if (found) {
                rule = 0;
            }
        }
    }
}
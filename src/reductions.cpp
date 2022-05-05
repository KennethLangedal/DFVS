#include "reductions.hpp"
#include <cassert>
#include <functional>

graph_search::graph_search(size_t N) : search(num_reductions), visited(num_reductions, bitvector(N)), DFS_visited(N), DFS_stack{}, SCC_id(N, N), L{}, tmp(N), tmp1(N), tmp2(N) {
    for (size_t r = 0; r < num_reductions; ++r) {
        for (size_t u = 0; u < N; ++u) {
            search[r].push_back(u);
        }
    }
}

void push_search(graph_search &gs, size_t u) {
    for (size_t r = 0; r < num_reductions; ++r) {
        if (gs.visited[r].get(u))
            gs.search[r].push_back(u);
        gs.visited[r].reset(u);
    }
}

size_t pop_search(graph_search &gs, size_t rule) {
    assert(rule < num_reductions && !gs.search[rule].empty());
    size_t u = gs.search[rule].back();
    gs.search[rule].pop_back();
    gs.visited[rule].set(u);
    return u;
}

bool zero_degree_reduction(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    assert(g.active_vertices().get(u));
    if (g.in_degree(u) == 0 || g.out_degree(u) == 0) {
        exclude_from_fvs(g, fvs, gs, u);
        return true;
    }
    return false;
}

bool self_edge_reduction(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    assert(g.active_vertices().get(u));
    if (g.out(u).get(u)) {
        add_to_fvs(g, fvs, gs, u);
        return true;
    }
    return false;
}

bool one_degree_reduction(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    assert(g.active_vertices().get(u));
    if (g.in_degree(u) == 1 || g.out_degree(u) == 1) {
        exclude_from_fvs(g, fvs, gs, u);
        return true;
    }
    return false;
}

bool redundant_edges(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    assert(g.active_vertices().get(u));

    // Zero in or out (excluding 2-cycles)

    gs.tmp.set_and(g.out(u), g.in(u));
    if (gs.tmp == g.out(u) && gs.tmp != g.in(u)) { // if ((g.out(u) - g.in(u)).empty() && !((g.in(u) - g.out(u)).empty()))
        gs.tmp.set_and_not(g.in(u), g.out(u));     // Make tmp copy for large version, could invalidate iterator
        for (size_t v : gs.tmp) {
            g.remove_edge(v, u);
            push_search(gs, v);
        }
        push_search(gs, u);
        return true;
    }
    if (gs.tmp == g.in(u) && gs.tmp != g.out(u)) {
        gs.tmp.set_and_not(g.out(u), g.in(u));
        for (size_t v : gs.tmp) {
            g.remove_edge(u, v);
            push_search(gs, v);
        }
        push_search(gs, u);
        return true;
    }

    // DOME rule

    bool res = false;
    gs.tmp.set_and_not(g.out(u), g.in(u));
    gs.tmp1.set_and_not(g.in(u), g.out(u));
    for (size_t v : gs.tmp) {
        gs.tmp2.set_and_not(g.out(v), g.in(v));
        if (gs.tmp1.subset_eq(g.in(v)) || gs.tmp2.subset_eq(g.out(u))) {
            g.remove_edge(u, v);
            push_search(gs, v);
            res = true;
        }
    }

    if (res)
        push_search(gs, u);

    return res;
}

bool isolated_vertex_reduction(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    assert(g.active_vertices().get(u));
    if (g.out(u) == g.in(u)) {
        for (size_t v : g.out(u)) {
            if (g.out(u).intersection_size(g.out(v)) != g.out_degree(u) - 1) // Not part of clique
                return false;
        }
        exclude_from_fvs(g, fvs, gs, u);
        return true;
    }
    return false;
}

bool twin_vertices_reduction(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    assert(g.active_vertices().get(u));

    if (g.out(u) == g.in(u) && !g.out(u).get(u)) {
        gs.tmp.clear();
        gs.tmp.set(u);
        for (size_t v : g.out(*g.out(u).begin())) {
            if (v != u && g.out(v) == g.in(v) && g.out(v) == g.out(u))
                gs.tmp.set(v);
        }
        if (gs.tmp.popcount() >= g.out_degree(u)) {
            for (size_t v : g.out(u))
                add_to_fvs(g, fvs, gs, v);
            for (size_t v : gs.tmp)
                exclude_from_fvs(g, fvs, gs, v);
            return true;
        } else if (gs.tmp.popcount() + 1 >= g.out_degree(u)) {
            for (size_t v : g.out(u)) {
                gs.tmp1.set_and(g.out(v), g.in(v));
                if (gs.tmp1.intersection_size(g.out(u)) > 0) {
                    for (size_t v : g.out(u))
                        add_to_fvs(g, fvs, gs, v);
                    for (size_t v : gs.tmp)
                        exclude_from_fvs(g, fvs, gs, v);
                    return true;
                }
            }
        }
    }
    return false;
}

bool dominating_vertex_reduction(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    assert(g.active_vertices().get(u));

    gs.tmp.set_and(g.out(u), g.in(u));
    for (size_t v : gs.tmp) {
        if (gs.tmp.intersection_size(g.out(v)) == g.out_degree(v) - 1 || gs.tmp.intersection_size(g.in(v)) == g.in_degree(v) - 1) {
            add_to_fvs(g, fvs, gs, u);
            return true;
        }
    }

    for (size_t v : gs.tmp) {
        gs.tmp1.set_and(g.out(v), g.in(v));
        if (gs.tmp1.intersection_size(g.out(u)) == g.out_degree(u) - 1 || gs.tmp1.intersection_size(g.in(u)) == g.in_degree(u) - 1) {
            add_to_fvs(g, fvs, gs, v);
            return true;
        }
    }
    return false;
}

bool cycle_dominating_vertex(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    assert(g.active_vertices().get(u));

    gs.tmp.set_and_not(g.out(u), g.in(u));
    gs.tmp2.set_and(g.out(u), g.in(u));
    if (gs.tmp.popcount() > 0) {
        for (size_t v : gs.tmp) {
            gs.tmp1.set_and_not(g.out(v), g.in(v));
            for (size_t w : gs.tmp1) {
                if (g.out(w).get(u) && !g.in(w).get(u)) { // length 3 cycle
                    if ((gs.tmp2.intersection_size(g.out(v)) == g.out_degree(v) - 1 && gs.tmp2.intersection_size(g.out(w)) == g.out_degree(w) - 1) ||
                        (gs.tmp2.intersection_size(g.in(v)) == g.in_degree(v) - 1 && gs.tmp2.intersection_size(g.in(w)) == g.in_degree(w) - 1)) {
                        add_to_fvs(g, fvs, gs, u);
                        return true;
                    }
                }
            }
        }
    }

    return false;
}

bool clique_and_one_fold(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    assert(g.active_vertices().get(u) && g.out_degree(u) >= 2);

    if (g.out(u) == g.in(u)) {
        size_t v = g.size();
        if (g.out_degree(u) == 2) {
            v = *g.out(u).begin();
            if (g.out(v) != g.in(v))
                return false;
        } else {
            for (size_t w : g.out(u)) {
                if (g.out(w) != g.in(w))
                    return false;
                gs.tmp.set_and(g.out(w), g.in(w));
                if (gs.tmp.intersection_size(g.out(u)) == 0) {
                    if (v != g.size())
                        return false;
                    v = w;
                } else if (gs.tmp.intersection_size(g.out(u)) != g.out_degree(u) - 2) {
                    return false;
                }
            }
            if (v == g.size())
                return false;
        }

        for (size_t w1 : g.out(u)) {
            push_search(gs, w1);
            for (size_t w2 : g.out(w1)) {
                push_search(gs, w2);
            }
            for (size_t w2 : g.in(w1)) {
                push_search(gs, w2);
            }
        }
        g.fold_clique_and_one(u, v);
        fvs.set(u);
        return true;
    }
    return false;
}

bool square_fold(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    assert(g.active_vertices().get(u));

    if (g.in(u) == g.out(u) && g.out_degree(u) == 4) {
        for (size_t v : g.out(u)) {
            gs.tmp.set_and(g.out(v), g.in(v));
            if (g.out(v) != g.in(v) || g.out(u).intersection_size(gs.tmp) != 2)
                return false;
        }

        size_t v1 = *g.out(u).begin(), v2;
        gs.tmp.set_and(g.out(v1), g.in(v1));
        gs.tmp1.set_and(gs.tmp, g.out(u));
        auto it = gs.tmp1.begin();
        size_t w1 = *it;
        ++it;
        size_t w2 = *it;
        for (size_t v : g.out(u)) {
            if (v != v1 && v != w1 && v != w2)
                v2 = v;
        }

        g.fold_square(u, v1, v2, w1, w2);
        fvs.set(v2);
        fvs.set(w2);
        for (size_t v : g.out(v1))
            push_search(gs, v);
        for (size_t v : g.in(v1))
            push_search(gs, v);
        for (size_t v : g.out(w1))
            push_search(gs, v);
        for (size_t v : g.in(w1))
            push_search(gs, v);
        return true;
    }
    return false;
}

bool has_cycle(const graph &g, bitvector &visited, size_t u, size_t &length) {
    static bitvector current(g.size()), next(g.size()), tmp(g.size());
    if (current.size() != g.size()) {
        current = bitvector(g.size()), next = bitvector(g.size()), tmp = bitvector(g.size());
    }
    current.set_and_not(g.out(u), g.in(u));
    current.set_and_not(current, visited);
    visited |= current;
    length = 0;
    while (current.popcount() > 0) {
        ++length;
        next.clear();
        if (current.get(u))
            return true;
        for (size_t v : current) {
            next |= tmp.set_and_not(g.out(v), g.in(v));
        }
        next.set_and_not(next, visited);
        visited |= next;
        std::swap(next, current);
    }
    return false;
}

bool specific_pattern(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    assert(g.active_vertices().get(u));

    gs.tmp.set_and(g.in(u), g.out(u));

    if (gs.tmp.popcount() == 0) {
        bool res = false;

        // 2 in and 2 out, any double edge between neighbours gives reduction
        if (g.out_degree(u) == 2 && g.in_degree(u) == 2) {
            gs.tmp1.set_or(g.out(u), g.in(u));
            for (size_t v : gs.tmp1) {
                gs.tmp2.set_and(g.out(v), g.in(v));
                if (gs.tmp2.intersection_size(gs.tmp1) > 0)
                    res = true;
            }
            // 2 out, double edge between them gives reduction
        } else if (g.out_degree(u) == 2) {
            size_t n[2], i = 0;
            for (size_t v : g.out(u))
                n[i++] = v;
            if (g.out(n[0]).get(n[1]) && g.out(n[1]).get(n[0]))
                res = true;
            // 2 in, double edge between them gives reduction
        } else if (g.in_degree(u) == 2) {
            size_t n[2], i = 0;
            for (size_t v : g.in(u))
                n[i++] = v;
            if (g.out(n[0]).get(n[1]) && g.out(n[1]).get(n[0]))
                res = true;
        }
        if (res) {
            exclude_from_fvs(g, fvs, gs, u);
            return true;
        }
    }

    // Single double edge, where the other also connects the same vertices

    if (gs.tmp.popcount() == 1) {
        size_t v = *gs.tmp.begin();
        gs.tmp1.set_and_not(g.out(u), g.in(u));
        gs.tmp2.set_and_not(g.in(u), g.out(u));
        if (gs.tmp1.subset_eq(g.out(v)) && gs.tmp2.subset_eq(g.in(v))) {
            exclude_from_fvs(g, fvs, gs, u);
            return true;
        }
    }

    // Cycle blocking tests

    size_t length;
    // 2 cycle neighbour
    if (gs.tmp.popcount() == 1) {
        size_t v = *gs.tmp.begin();
        gs.tmp1.set_and(g.out(v), g.in(v));
        if (gs.tmp1.popcount() > 1) {
            gs.tmp1.reset(u);
            if (!has_cycle(g, gs.tmp1, u, length)) {
                add_to_fvs(g, fvs, gs, v);
                return true;
            }
        }
    }

    return false;
}

bool SCC_edge_reduction(graph &g, bitvector &fvs, const bitvector &nodes, graph_search &gs, std::vector<bitvector> &SCC) {
    gs.DFS_stack.clear();
    gs.L.clear();
    gs.DFS_visited.clear();
    SCC.clear();

    std::function<void(size_t)> dfs_visit = [&](size_t u) {
        if (gs.DFS_visited.get(u))
            return;
        gs.DFS_visited.set(u);
        for (size_t v : g.in(u))
            dfs_visit(v);
        gs.L.push_back(u);
    };

    for (size_t u : nodes) {
        gs.SCC_id[u] = g.size();
        dfs_visit(u);
    }

    std::function<void(size_t, size_t, bitvector &)> assign = [&](size_t u, size_t root, bitvector &CC) {
        if (gs.SCC_id[u] != g.size())
            return;
        CC.set(u);
        gs.SCC_id[u] = root;
        for (size_t v : g.out(u))
            assign(v, root, CC);
    };

    for (size_t u : gs.L) {
        if (gs.SCC_id[u] != g.size())
            continue;
        bitvector CC(g.size());
        assign(u, u, CC);
        SCC.push_back(CC);
    }

    bool res = false;
    for (size_t u : nodes) {
        for (size_t v : g.out(u)) {
            if (gs.SCC_id[u] != gs.SCC_id[v]) {
                g.remove_edge(u, v);
                push_search(gs, u);
                push_search(gs, v);
                res = true;
            }
        }
    }

    return res;
}

bool global_redundant_edge(graph &g, bitvector &fvs, graph_search &gs) {
    // gs.tmp = path
    std::function<bool(size_t, size_t, size_t)> dfs = [&](size_t u, size_t v, size_t d) {
        if (gs.tmp.intersection_size(g.in(v)) > 1)
            return false;
        if (g.out(v).get(u)) {
            return gs.tmp.intersection_size(g.out(v)) == 1;
        }
        if (gs.tmp.intersection_size(g.out(v)) > 0) {
            return false;
        }
        if (d > 20)
            return true;
        bool res = false;
        gs.tmp.set(v);
        bool single = g.out(u).popcount() - g.out(u).intersection_size(g.in(u)) == 1;
        for (size_t w : g.out(v)) {
            res |= dfs(u, w, d + (single ? 0 : 1));
            if (res)
                break;
        }
        gs.tmp.reset(v);
        return res;
    };

    gs.tmp.clear();

    bool res = false;
    for (size_t u : g.active_vertices()) {
        gs.tmp.set(u);
        for (size_t v : g.out(u)) {
            if (!dfs(u, v, 0)) {
                g.remove_edge(u, v);
                push_search(gs, u);
                push_search(gs, v);
                res = true;
            }
        }
        gs.tmp.reset(u);
    }
    return res;
}

void add_to_fvs(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    fvs.set(u);
    for (size_t v : g.out(u))
        push_search(gs, v);
    for (size_t v : g.in(u))
        push_search(gs, v);
    g.remove_vertex(u);
}

void exclude_from_fvs(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    for (size_t v : g.out(u))
        push_search(gs, v);
    for (size_t v : g.in(u))
        push_search(gs, v);
    gs.tmp.set_and_not(g.in(u), g.out(u));
    gs.tmp1.set_and_not(g.out(u), g.in(u));
    gs.tmp2.set_and(g.out(u), g.in(u));

    for (size_t v : gs.tmp2)
        g.add_edge(v, v);
    for (size_t v1 : gs.tmp) {
        for (size_t v2 : gs.tmp1) {
            g.add_edge(v1, v2);
        }
    }
    g.remove_vertex(u);
}

void reduce_graph(graph &g, bitvector &fvs, const bitvector &nodes, graph_search &gs, std::vector<bitvector> &SCC, bool global) {
    size_t rule = 0;
    while (rule < num_reductions) {
        if (gs.search[rule].empty()) {
            rule++;
        } else {
            size_t u = pop_search(gs, rule);
            if (!g.active_vertices().get(u))
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
            case reductions::cycle_dominating_vertex:
                found = cycle_dominating_vertex(g, fvs, gs, u);
                break;
            case reductions::twin_vertices:
                found = twin_vertices_reduction(g, fvs, gs, u);
                break;
            case reductions::clique_and_one_fold:
                found = clique_and_one_fold(g, fvs, gs, u);
                break;
            case reductions::square_fold:
                found = square_fold(g, fvs, gs, u);
                break;
            case reductions::specific_pattern:
                found = specific_pattern(g, fvs, gs, u);
                break;
            default:
                break;
            }
            if (found) {
                rule = 0;
            }
        }

        if (rule == num_reductions) {
            gs.tmp2.set_and(nodes, g.active_vertices());
            if (SCC_edge_reduction(g, fvs, gs.tmp2, gs, SCC)) {
                rule = 0;
            } else if ((global || g.active_vertices().popcount() < 100) && global_redundant_edge(g, fvs, gs)) {
                rule = 0;
            }
        }
    }
}
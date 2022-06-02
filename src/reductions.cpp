#include "reductions.hpp"
#include "flow_graph.hpp"
#include <algorithm>
#include <cassert>

graph_search::graph_search(size_t N, bool queue_all)
    : search(num_reductions), visited(num_reductions, bitvector(N)), SCC{}, DFS_visited(N), DFS_stack{}, SCC_id(N, N), L{}, tmp{} {
    if (!queue_all) {
        for (size_t r = 0; r < num_reductions; ++r) {
            visited[r].fill();
        }
        return;
    }
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

void queue_neighbourhood(sparse_graph &g, graph_search &gs, size_t u) {
    for (auto v : g.out(u)) {
        push_search(gs, v);
        for (auto w : g.out(v))
            push_search(gs, w);
        for (auto w : g.in(v))
            push_search(gs, w);
    }
    for (auto v : g.in(u)) {
        push_search(gs, v);
        for (auto w : g.out(v))
            push_search(gs, w);
        for (auto w : g.in(v))
            push_search(gs, w);
    }
}

void add_to_fvs(sparse_graph &g, bitvector &fvs, graph_search &gs, reduction_engine &re, size_t u) {
    fvs.set(u);
    queue_neighbourhood(g, gs, u);
    re.remove_include_vertex(g, u);
}

void exclude_from_fvs(sparse_graph &g, bitvector &fvs, graph_search &gs, reduction_engine &re, size_t u) {
    queue_neighbourhood(g, gs, u);
    re.remove_exclude_vertex(g, u);
}

// ************
// *Reductions*
// ************

bool zero_degree_reduction(sparse_graph &g, bitvector &fvs, graph_search &gs, reduction_engine &re, size_t u) {
    assert(g.is_active(u));
    if (g.in_degree(u) == 0 || g.out_degree(u) == 0) {
        exclude_from_fvs(g, fvs, gs, re, u);
        return true;
    }
    return false;
}

bool self_edge_reduction(sparse_graph &g, bitvector &fvs, graph_search &gs, reduction_engine &re, size_t u) {
    assert(g.is_active(u));
    if (g.has_edge(u, u)) {
        add_to_fvs(g, fvs, gs, re, u);
        return true;
    }
    return false;
}

bool one_degree_reduction(sparse_graph &g, bitvector &fvs, graph_search &gs, reduction_engine &re, size_t u) {
    assert(g.is_active(u));
    if (g.in_degree(u) == 1 || g.out_degree(u) == 1) {
        exclude_from_fvs(g, fvs, gs, re, u);
        return true;
    }
    return false;
}

bool redundant_edges(sparse_graph &g, bitvector &fvs, graph_search &gs, reduction_engine &re, size_t u) {
    assert(g.is_active(u));

    // Zero in or out (excluding 2-cycles)

    if (g.out_degree_non_pi(u) == 0 && g.in_degree_non_pi(u) > 0) {
        gs.tmp.clear();
        for (auto v : g.in_non_pi(u))
            gs.tmp.push_back(v);

        for (auto v : gs.tmp) {
            re.remove_edge(g, v, u);
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
            re.remove_edge(g, u, v);
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
            re.remove_edge(g, u, v);
            push_search(gs, v);
        }
        push_search(gs, u);
    }

    return res;
}

bool isolated_vertex_reduction(sparse_graph &g, bitvector &fvs, graph_search &gs, reduction_engine &re, size_t u) {
    assert(g.is_active(u));
    if (g.pi_degree(u) == g.out_degree(u)) {
        for (auto v : g.pi(u)) {
            if (intersection_size(g.out(u), g.out(v)) != g.out_degree(u) - 1) // Not part of clique
                return false;
        }
        exclude_from_fvs(g, fvs, gs, re, u);
        return true;
    }
    return false;
}

bool dominating_vertex_reduction(sparse_graph &g, bitvector &fvs, graph_search &gs, reduction_engine &re, size_t u) {
    assert(g.is_active(u));

    for (auto v : g.pi(u)) {
        if (intersection_size(g.pi(u), g.out(v)) == g.out_degree(v) - 1 || intersection_size(g.pi(u), g.in(v)) == g.in_degree(v) - 1) {
            add_to_fvs(g, fvs, gs, re, u);
            return true;
        }
        if (intersection_size(g.pi(v), g.out(u)) == g.out_degree(u) - 1 || intersection_size(g.pi(v), g.in(u)) == g.in_degree(u) - 1) {
            add_to_fvs(g, fvs, gs, re, v);
            return true;
        }
    }

    // unconfined rule, VC only

    if (g.pi_only(u)) {
        static std::vector<uint32_t> S, NS, NSI;
        S.clear();
        NS.clear();
        NSI.clear();
        S.push_back(u);
        std::copy(std::begin(g.pi(u)), std::end(g.pi(u)), std::back_inserter(NS));
        std::copy(std::begin(g.pi(u)), std::end(g.pi(u)), std::back_inserter(NSI));
        NSI.insert(std::lower_bound(std::begin(NSI), std::end(NSI), u), u);
        while (true) {
            uint32_t w = g.size();
            for (auto v : NS) {
                if (g.pi_only(v) && intersection_size(g.pi(v), S) == 1) {
                    gs.tmp.clear();
                    std::set_difference(std::begin(g.pi(v)), std::end(g.pi(v)), std::begin(NSI), std::end(NSI), std::back_inserter(gs.tmp));
                    if (gs.tmp.empty()) {
                        add_to_fvs(g, fvs, gs, re, u);
                        return true;
                    } else if (gs.tmp.size() == 1 && g.pi_only(gs.tmp.front())) {
                        w = gs.tmp.front();
                    }
                }
            }
            if (w == g.size())
                break;
            S.insert(std::lower_bound(std::begin(S), std::end(S), w), w);
            gs.tmp.clear();
            std::set_union(std::begin(NS), std::end(NS), std::begin(g.pi(w)), std::end(g.pi(w)), std::back_inserter(gs.tmp));
            std::swap(gs.tmp, NS);

            NSI.insert(std::lower_bound(std::begin(NSI), std::end(NSI), w), w);
            gs.tmp.clear();
            std::set_union(std::begin(NSI), std::end(NSI), std::begin(g.pi(w)), std::end(g.pi(w)), std::back_inserter(gs.tmp));
            std::swap(gs.tmp, NSI);
        }
    }

    return false;
}

bool cycle_dominating_vertex(sparse_graph &g, bitvector &fvs, graph_search &gs, reduction_engine &re, size_t u) {
    assert(g.is_active(u));

    if (g.out_degree_non_pi(u) > 10)
        return false;

    for (auto v : g.out_non_pi(u)) {
        if (g.out_degree_non_pi(v) > 10)
            continue;
        for (auto w : g.out_non_pi(v)) {
            if (g.has_edge(w, u) && !g.has_edge(u, w)) { // length 3 cycle
                if ((intersection_size(g.pi(u), g.out(v)) == g.out_degree(v) - 1 && intersection_size(g.pi(u), g.out(w)) == g.out_degree(w) - 1) ||
                    (intersection_size(g.pi(u), g.in(v)) == g.in_degree(v) - 1 && intersection_size(g.pi(u), g.in(w)) == g.in_degree(w) - 1)) {
                    add_to_fvs(g, fvs, gs, re, u);
                    return true;
                }
            }
        }
    }

    return false;
}

bool twin_vertices_reduction(sparse_graph &g, bitvector &fvs, graph_search &gs, reduction_engine &re, size_t u) {
    assert(g.is_active(u) && g.out_degree(u) > 0);

    gs.tmp.clear();
    gs.tmp.push_back(u);

    for (auto v : g.in(*g.out(u).begin())) {
        if (v != u && g.out(u) == g.out(v) && g.in(u) == g.in(v)) {
            gs.tmp.push_back(v);
        }
    }

    if (gs.tmp.size() >= g.out_degree(u) || gs.tmp.size() >= g.in_degree(u)) {
        for (auto v : gs.tmp)
            exclude_from_fvs(g, fvs, gs, re, v);
        return true;
    }

    if (g.pi_only(u) && g.pi_degree(u) == 3 && gs.tmp.size() == 2) {
        for (auto v : g.pi(u)) {
            if (intersection_size(g.pi(u), g.pi(v)) > 0) {
                for (auto v : gs.tmp)
                    exclude_from_fvs(g, fvs, gs, re, v);
                return true;
            }
        }
        for (auto v : g.pi(u)) {
            if (!g.pi_only(v))
                return false;
        }
        uint32_t f = *std::begin(g.pi(u));
        re.fold_twin(g, gs.tmp[0], gs.tmp[1], fvs);
        queue_neighbourhood(g, gs, f);
        push_search(gs, f);
        return true;
    }

    return false;
}

bool clique_and_one_fold(sparse_graph &g, bitvector &fvs, graph_search &gs, reduction_engine &re, size_t u) {
    assert(g.is_active(u) && g.out_degree(u) >= 2);

    if (g.out_degree_non_pi(u) == 0 && g.in_degree_non_pi(u) == 0) {
        auto v = g.size();
        if (g.out_degree(u) == 2) {
            for (auto w : g.out(u)) {
                if (g.out_degree_non_pi(u) == 0 && g.in_degree_non_pi(u) == 0)
                    v = w;
            }
        } else {
            for (auto w : g.out(u)) {
                if (intersection_size(g.pi(w), g.out(u)) == 0) {
                    if (v != g.size())
                        return false;
                    v = w;
                } else if (intersection_size(g.pi(w), g.out(u)) != g.out_degree(u) - 2) {
                    return false;
                }
            }
        }
        if (v == g.size() || g.out_degree_non_pi(v) > 0 || g.in_degree_non_pi(v) > 0)
            return false;

        for (auto w1 : g.out(u)) {
            push_search(gs, w1);
            queue_neighbourhood(g, gs, w1);
        }
        re.fold_clique_and_one(g, u, v, fvs);
        return true;
    }
    return false;
}

bool square_fold(sparse_graph &g, bitvector &fvs, graph_search &gs, reduction_engine &re, size_t u) {
    assert(g.is_active(u));

    if (g.out_degree_non_pi(u) == 0 && g.in_degree_non_pi(u) == 0 && g.pi_degree(u) == 4) {
        for (auto v : g.pi(u)) {
            if (g.out_degree_non_pi(v) != 0 || g.in_degree_non_pi(v) == 0 || intersection_size(g.pi(u), g.pi(v)) != 2)
                return false;
        }

        auto it = g.pi(u).begin();
        uint32_t v1 = *it++, v2, w1, w2;
        if (g.has_edge(v1, *it)) {
            w1 = *it++;
            if (g.has_edge(v1, *it)) {
                w2 = *it++;
                v2 = *it;
            } else {
                v2 = *it++;
                w2 = *it;
            }
        } else {
            v2 = *it++;
            w1 = *it++;
            w2 = *it;
        }

        re.fold_square(g, u, v1, v2, w1, w2, fvs);
        queue_neighbourhood(g, gs, v1);
        queue_neighbourhood(g, gs, w1);
        return true;
    }
    return false;
}

bool specific_pattern(sparse_graph &g, bitvector &fvs, graph_search &gs, reduction_engine &re, size_t u) {
    assert(g.is_active(u));

    if (g.pi_degree(u) == 0) {
        bool res = false;

        // 2 in and 2 out, any double edge between neighbours gives reduction
        if (g.out_degree(u) == 2 && g.in_degree(u) == 2) {
            for (auto v : g.in(u)) {
                for (auto w : g.out(u)) {
                    if (g.has_edge(v, w) && g.has_edge(w, v))
                        res = true;
                }
            }
        }
        // 2 out, double edge between them gives reduction
        if (g.out_degree(u) == 2) {
            size_t n[2], i = 0;
            for (size_t v : g.out(u))
                n[i++] = v;
            if (g.has_edge(n[0], n[1]) && g.has_edge(n[1], n[0]))
                res = true;
        }
        // 2 in, double edge between them gives reduction
        if (g.in_degree(u) == 2) {
            size_t n[2], i = 0;
            for (size_t v : g.in(u))
                n[i++] = v;
            if (g.has_edge(n[0], n[1]) && g.has_edge(n[1], n[0]))
                res = true;
        }
        if (res) {
            exclude_from_fvs(g, fvs, gs, re, u);
            return true;
        }
    }

    // Single double edge, where the other also connects the same vertices

    if (g.pi_degree(u) == 1) {
        size_t v = *g.pi(u).begin();
        if (includes(g.out(v), g.out_non_pi(u)) && includes(g.in(v), g.in_non_pi(u))) {
            exclude_from_fvs(g, fvs, gs, re, u);
            return true;
        }
    }

    return false;
}

// ************
// *global-SCC*
// ************

bool SCC_edge_reduction(sparse_graph &g, bitvector &fvs, graph_search &gs, reduction_engine &re) {
    gs.DFS_stack.clear();
    gs.L.clear();
    gs.DFS_visited.clear();
    gs.SCC.clear();

    std::function<void(uint32_t)> dfs_visit = [&](auto u) {
        if (gs.DFS_visited.get(u))
            return;
        gs.DFS_visited.set(u);
        for (auto v : g.in(u))
            dfs_visit(v);
        gs.L.push_back(u);
    };

    for (auto u : g.active_vertices()) {
        gs.SCC_id[u] = g.size();
        dfs_visit(u);
    }

    std::function<void(uint32_t, uint32_t, bitvector &)> assign = [&](auto u, auto root, bitvector &CC) {
        if (gs.SCC_id[u] != g.size())
            return;
        CC.set(u);
        gs.SCC_id[u] = root;
        for (auto v : g.out(u))
            assign(v, root, CC);
    };

    for (auto u : gs.L) {
        if (gs.SCC_id[u] != g.size())
            continue;
        bitvector CC(g.size());
        assign(u, u, CC);
        gs.SCC.push_back(CC);
    }

    bool res = false;
    for (auto u : g.active_vertices()) {
        for (auto v : g.out(u)) {
            if (gs.SCC_id[u] != gs.SCC_id[v]) {
                re.remove_edge(g, u, v);
                push_search(gs, u);
                push_search(gs, v);
                res = true;
            }
        }
    }

    return res;
}

bool reduction_half_lp(sparse_graph &g, bitvector &fvs, graph_search &gs, reduction_engine &re) {
    for (auto u : g.active_vertices()) {
        if (!g.pi_only(u))
            return false;
    }

    std::vector<std::tuple<int32_t, int32_t, int32_t>> edges;
    int32_t s = g.size() * 2, t = (g.size() * 2) + 1;
    for (auto u : g.active_vertices()) {
        edges.push_back({s, u, 1});
        edges.push_back({g.size() + u, t, 1});
        for (auto v : g.out(u)) {
            edges.push_back({u, g.size() + v, 1});
        }
    }

    flow_graph<int32_t, int32_t> fg(t + 1, edges);
    int32_t flow = fg.solve(s, t), active_count = g.active_vertices().popcount();
    std::vector<uint32_t> current, next;
    std::vector<bool> visited(g.size() * 2, false), L(g.size(), true), R(g.size(), false);
    for (auto &&[u, w] : fg[s]) {
        if (w > 0) {
            current.push_back(u);
            visited[u] = true;
        }
    }

    while (current.size() > 0) {
        next.clear();
        for (auto u : current) {
            if (u < g.size()) {
                L[u] = false;
            } else {
                R[u - g.size()] = true;
            }
            for (auto [v, w] : fg[u]) {
                if (v < g.size() * 2 && w > 0 && !visited[v]) {
                    next.push_back(v);
                    visited[v] = true;
                }
            }
        }
        std::swap(current, next);
    }

    bool res = false;
    for (auto u : g.active_vertices()) {
        if (g.is_active(u) && !L[u] && !R[u]) {
            res = true;
            exclude_from_fvs(g, fvs, gs, re, u);
        } else if (g.is_active(u) && L[u] && R[u]) {
            res = true;
            add_to_fvs(g, fvs, gs, re, u);
        }
    }

#ifdef VERBOSE
    if (res)
        std::cout << "LP did something" << std::endl;
#endif

    return res;
}

void reduce_graph(sparse_graph &g, bitvector &fvs, graph_search &gs, reduction_engine &re, bool SCC) {
    size_t rule = 0;
    bool found;
    while (rule < num_reductions) {
        found = false;
        if (gs.search[rule].empty()) {
            rule++;
        } else {
            size_t u = pop_search(gs, rule);
            if (!g.is_active(u))
                continue;
            switch ((reductions)rule) {
            case reductions::zero_degree:
                found = zero_degree_reduction(g, fvs, gs, re, u);
                break;
            case reductions::self_edge:
                found = self_edge_reduction(g, fvs, gs, re, u);
                break;
            case reductions::one_degree:
                found = one_degree_reduction(g, fvs, gs, re, u);
                break;
            case reductions::redundant_edges:
                found = redundant_edges(g, fvs, gs, re, u);
                break;
            case reductions::isolated_vertex:
                found = isolated_vertex_reduction(g, fvs, gs, re, u);
                break;
            case reductions::dominating_vertex:
                found = dominating_vertex_reduction(g, fvs, gs, re, u);
                break;
            case reductions::cycle_dominating_vertex:
                found = cycle_dominating_vertex(g, fvs, gs, re, u);
                break;
            case reductions::twin_vertices:
                found = twin_vertices_reduction(g, fvs, gs, re, u);
                break;
            case reductions::clique_and_one_fold:
                found = clique_and_one_fold(g, fvs, gs, re, u);
                break;
            case reductions::square_fold:
                found = square_fold(g, fvs, gs, re, u);
                break;
            case reductions::specific_pattern:
                found = specific_pattern(g, fvs, gs, re, u);
                break;
            default:
                break;
            }
        }
        if (found || (SCC && rule == num_reductions && (SCC_edge_reduction(g, fvs, gs, re) || reduction_half_lp(g, fvs, gs, re)))) {
            rule = 0;
        }
    }
}
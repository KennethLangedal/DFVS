#include "reductions.hpp"

graph_search::graph_search() : search(num_reductions), visited(num_reductions, bitvector<N>{}), DFS_visited{}, DFS_stack{}, SCC_id(N * 64, N * 64), L{} {
    for (size_t r = 0; r < num_reductions; ++r) {
        for (size_t u = 0; u < N * 64; ++u) {
            search[r].push_back(u);
        }
    }
}

void push_search(graph_search &gs, size_t u) {
    assert(u < N * 64);
    for (size_t r = 0; r < num_reductions; ++r) {
        if (test(gs.visited[r], u))
            gs.search[r].push_back(u);
        reset(gs.visited[r], u);
    }
}

size_t pop_search(graph_search &gs, size_t rule) {
    assert(rule < num_reductions && !gs.search[rule].empty());
    size_t u = gs.search[rule].back();
    gs.search[rule].pop_back();
    set(gs.visited[rule], u);
    return u;
}

bool zero_degree_reduction(graph &g, bitvector<N> &fvs, graph_search &gs, size_t u) {
    assert(u < N * 64 && test(g.active_vertices(), u));
    if (g.in_degree(u) == 0 || g.degree(u) == 0) {
        deactivate_vertex(g, gs, u);
        return true;
    }
    return false;
}

bool self_edge_reduction(graph &g, bitvector<N> &fvs, graph_search &gs, size_t u) {
    assert(u < N * 64 && test(g.active_vertices(), u));
    if (g.self_loop(u)) {
        add_to_fvs(g, fvs, gs, u);
        return true;
    }
    return false;
}

bool one_degree_reduction(graph &g, bitvector<N> &fvs, graph_search &gs, size_t u) {
    assert(u < N * 64 && test(g.active_vertices(), u));
    if (g.in_degree(u) == 1 || g.degree(u) == 1) {
        exclude_from_fvs(g, fvs, gs, u);
        return true;
    }
    return false;
}

bool redundant_edges(graph &g, bitvector<N> &fvs, graph_search &gs, size_t u) {
    assert(u < N * 64 && test(g.active_vertices(), u));
    bitvector<N> tmp = g.out(u) & g.in(u);
    if (tmp == g.out(u) && tmp != g.in(u)) {
        visit(g.in(u) & ~tmp, [&](size_t v) {
            g.remove_edge(v, u);
            push_search(gs, v);
        });
        push_search(gs, u);
        return true;
    }
    if (tmp == g.in(u) && tmp != g.out(u)) {
        visit(g.out(u) & ~tmp, [&](size_t v) {
            g.remove_edge(u, v);
            push_search(gs, v);
        });
        push_search(gs, u);
        return true;
    }

    bool res = false;
    bitvector<N> single_out = g.out(u) & ~g.in(u), single_in = g.in(u) & ~g.out(u);
    size_t count = popcount(single_out);
    visit(single_in, [&](size_t v) {
        if (intersection_size(g.out(v), single_out) == count) {
            g.remove_edge(v, u);
            push_search(gs, v);
            res = true;
        }
    });
    if (res)
        push_search(gs, u);

    return res;
}

bool isolated_vertex_reduction(graph &g, bitvector<N> &fvs, graph_search &gs, size_t u) {
    assert(u < N * 64 && test(g.active_vertices(), u));
    bitvector<N> tmp = g.out(u) & g.in(u);
    if (tmp == g.out(u) || tmp == g.in(u)) {
        bool res = true;
        size_t degree = popcount(tmp);
        visit(tmp, [&](size_t v) { res &= intersection_size(tmp, g.out(v)) == degree - 1; });
        if (res) {
            visit(tmp, [&](size_t v) { add_to_fvs(g, fvs, gs, v); });
            deactivate_vertex(g, gs, u);
            return true;
        }
    }
    return false;
}

bool twin_vertices_reduction(graph &g, bitvector<N> &fvs, graph_search &gs, size_t u) {
    assert(u < N * 64 && test(g.active_vertices(), u));
    bitvector<N> tmp = g.out(u) & g.in(u);
    if (tmp == g.out(u) || tmp == g.in(u)) {
        bitvector<N> twins{};
        set(twins, u);
        visit(g.out(first(tmp)), [&](size_t v) {
            if (v == u)
                return;
            bitvector<N> tmp2 = g.out(v) & g.in(v);
            if (tmp2 == tmp && (tmp2 == g.out(v) || tmp2 == g.in(v))) {
                set(twins, v);
            }
        });
        if (popcount(twins) >= popcount(tmp)) {
            visit(tmp, [&](size_t v) { add_to_fvs(g, fvs, gs, v); });
            visit(twins, [&](size_t v) { deactivate_vertex(g, gs, v); });
            return true;
        }
    }
    return false;
}

bool dominating_vertex_reduction(graph &g, bitvector<N> &fvs, graph_search &gs, size_t u) {
    assert(u < N * 64 && test(g.active_vertices(), u));
    bitvector<N> tmp = g.out(u) & g.in(u);
    if (tmp == g.out(u) || tmp == g.in(u)) {
        size_t degree = popcount(tmp);
        bool res = false;
        visit(tmp, [&](size_t v) {
            if (res)
                return;
            bitvector<N> tmp2 = g.out(v) & g.in(v);
            if (intersection_size(tmp, tmp2) == degree - 1) {
                add_to_fvs(g, fvs, gs, v);
                res = true;
            }
        });
        return res;
    }
    return false;
}

bool neighborhood_fold(graph &g, bitvector<N> &fvs, graph_search &gs, size_t u) {
    assert(u < N * 64 && test(g.active_vertices(), u));
    if (g.out(u) == g.in(u)) {
        if (popcount(g.out(u)) == 2) {
            g.fold_neighborhood(u);
            push_search(gs, u);
            visit(g.out(u), [&](size_t v) { push_search(gs, v); });
            visit(g.in(u), [&](size_t v) { push_search(gs, v); });
            return true;
        }
    }

    return false;
}

bool two_one_fold(graph &g, bitvector<N> &fvs, graph_search &gs, size_t u) {
    assert(u < N * 64 && test(g.active_vertices(), u));
    if (g.out(u) == g.in(u) && popcount(g.in(u)) == 3) {
        size_t v1 = N * 64, v2 = N * 64, w = N * 64;
        visit(g.out(u), [&](size_t x) {
            bitvector<N> in_and_out = g.out(x) & g.out(u);
            if (popcount(in_and_out) == 1) { // TODO, FIX BUG, can only have one doulbe edge and no single edges
                v1 = x;
                v2 = first(in_and_out);
            } else {
                w = x;
            }
        });
        if (v1 != N * 64 && v2 != N * 64 && w != N * 64) {
            visit(g.out(w), [&](size_t x) { push_search(gs, x); });
            visit(g.in(w), [&](size_t x) { push_search(gs, x); });
            g.fold_two_one(u, v1, v2, w);
            return true;
        }
    }

    return false;
}

bool specific_pattern(graph &g, bitvector<N> &fvs, graph_search &gs, size_t u) {
    assert(u < N * 64 && test(g.active_vertices(), u));

    bitvector<N> in_and_out = g.out(u) & g.in(u), out = g.out(u) & ~g.in(u), in = g.in(u) & ~g.out(u);
    if (popcount(in_and_out) == 1) {
        size_t v = first(in_and_out);
        bool covers_all_in = true, covers_all_out = true;
        visit(in, [&](size_t w) { covers_all_in = covers_all_in && test(g.in(v), w) && test(g.out(v), w); });
        visit(out, [&](size_t w) { covers_all_out = covers_all_out && test(g.in(v), w) && test(g.out(v), w); });

        if (covers_all_in || covers_all_out) {
            add_to_fvs(g, fvs, gs, v);
            return true;
        }
    }

    return false;
}

bool SCC_edge_reduction(graph &g, bitvector<N> &fvs, const bitvector<N> &nodes, graph_search &gs, std::vector<bitvector<N>> &SCC) {
    gs.DFS_stack.clear();
    gs.L.clear();
    gs.DFS_visited.fill(0ull);
    SCC.clear();

    std::function<void(size_t)> dfs_visit = [&](size_t u) {
        if (test(gs.DFS_visited, u))
            return;
        set(gs.DFS_visited, u);
        visit(g.out(u), [&](size_t v) { dfs_visit(v); });
        gs.L.push_back(u);
    };

    visit(nodes, [&](size_t u) {
        gs.SCC_id[u] = N * 64;
        dfs_visit(u);
    });

    std::function<void(size_t, size_t, bitvector<N> &)> assign = [&](size_t u, size_t root, bitvector<N> &CC) {
        if (gs.SCC_id[u] != N * 64)
            return;
        set(CC, u);
        gs.SCC_id[u] = root;
        visit(g.in(u), [&](size_t v) { assign(v, root, CC); });
    };

    for (auto u : gs.L) {
        if (gs.SCC_id[u] != N * 64)
            continue;
        bitvector<N> CC{};
        assign(u, u, CC);
        SCC.push_back(CC);
    }

    bool res = false;
    visit(nodes, [&](size_t u) {
        visit(g.out(u), [&](size_t v) {
            if (gs.SCC_id[u] != gs.SCC_id[v]) {
                g.remove_edge(u, v);
                push_search(gs, u);
                push_search(gs, v);
                res = true;
            }
        });
    });

    return res;
}

void add_to_fvs(graph &g, bitvector<N> &fvs, graph_search &gs, size_t u) {
    set(fvs, u);
    deactivate_vertex(g, gs, u);
}

void exclude_from_fvs(graph &g, bitvector<N> &fvs, graph_search &gs, size_t u) {
    bitvector<N> in = g.in(u) & ~g.out(u), out = g.out(u) & ~g.in(u), both = g.in(u)&g.out(u);
    visit(both, [&](size_t v) { g.add_edge(v, v); });
    visit(in, [&](size_t v1) { visit(out, [&](size_t v2) { g.add_edge(v1, v2); }); });
    deactivate_vertex(g, gs, u);
}

void deactivate_vertex(graph &g, graph_search &gs, size_t u) {
    visit(g.out(u), [&](size_t v) { push_search(gs, v); });
    visit(g.in(u), [&](size_t v) { push_search(gs, v); });
    g.remove_vertex(u);
}

void reduce_graph(graph &g, bitvector<N> &fvs, const bitvector<N> &nodes, graph_search &gs, std::vector<bitvector<N>> &SCC) {
    size_t rule = 0;
    while (rule < num_reductions) {
        if (gs.search[rule].empty()) {
            rule++;
        } else {
            size_t u = pop_search(gs, rule);
            if (!test(g.active_vertices(), u))
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
            case reductions::twin_vertices:
                found = twin_vertices_reduction(g, fvs, gs, u);
                break;
            case reductions::dominating_vertex:
                found = dominating_vertex_reduction(g, fvs, gs, u);
                break;
            case reductions::neighborhood_fold:
                found = neighborhood_fold(g, fvs, gs, u);
                break;
            case reductions::two_one_fold:
                found = two_one_fold(g, fvs, gs, u);
                break;
            case reductions::specific_pattern:
                found = specific_pattern(g, fvs, gs, u);
                break;
            default:
                break;
            }
            if (found)
                rule = 0;
        }

        if (rule == num_reductions && SCC_edge_reduction(g, fvs, nodes & g.active_vertices(), gs, SCC)) {
            rule = 0;
        }
    }
}
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
    if (g.self_loop(u)) {
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

    // Zero out or in (excluding 2-cycles)

    gs.tmp.set_and(g.out(u), g.in(u));
    if (gs.tmp == g.out(u) && gs.tmp != g.in(u)) {
        gs.tmp.set_and_not(g.in(u), g.out(u));
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
    if (g.out(u) == g.in(u)) {
        gs.tmp.clear();
        gs.tmp.set(u);
        for (size_t v : g.out(*g.out(u).begin())) {
            if (v != u && g.out(v) == g.in(v) && g.out(v) == g.in(u))
                gs.tmp.set(v);
        }
        if (gs.tmp.popcount() >= g.out_degree(u)) {
            for (size_t v : g.out(u))
                add_to_fvs(g, fvs, gs, v);
            for (size_t v : gs.tmp)
                exclude_from_fvs(g, fvs, gs, v);
            return true;
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

bool clique_and_one(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    assert(g.active_vertices().get(u) && g.out_degree(u) >= 2);
    if (g.out(u) == g.in(u)) {
        size_t v = g.size();
        if (g.out_degree(u) == 2) {
            v = *g.out(u).begin();
        } else {
            for (size_t w : g.out(u)) {
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
            for (size_t w2 : g.out(w1)) {
                push_search(gs, w2);
            }
        }
        g.fold_clique_and_one(u, v);
        return true;
    }
    return false;
}

bool square_fold(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    assert(g.active_vertices().get(u));

    if (g.in(u) == g.out(u) && g.out_degree(u) == 4) {
        for (size_t v : g.out(u)) {
            gs.tmp.set_and(g.out(v), g.in(v));
            if (g.out(u).intersection_size(gs.tmp) != 2)
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

// bool specific_pattern(graph &g, bitvector<N> &fvs, graph_search &gs, size_t u) {
//     assert(u < N * 64 && test(g.active_vertices(), u));

//     bitvector<N> in_and_out = g.out(u) & g.in(u), out = g.out(u) & ~g.in(u), in = g.in(u) & ~g.out(u);

//     // 2-cycle dominated
//     if (popcount(in_and_out) == 1) {
//         size_t v = first(in_and_out);
//         bool covers_all_in = true, covers_all_out = true;
//         visit(in, [&](size_t w) { covers_all_in = covers_all_in && test(g.in(v), w) && test(g.out(v), w); });
//         visit(out, [&](size_t w) { covers_all_out = covers_all_out && test(g.in(v), w) && test(g.out(v), w); });

//         if (covers_all_in || covers_all_out) {
//             add_to_fvs(g, fvs, gs, v);
//             return true;
//         }
//     }

//     if (popcount(in_and_out) == 0) {
//         bool res = false;

//         // 2 in and 2 out, any double edge between neighbours gives reduction
//         if (g.degree(u) == 2 && g.in_degree(u) == 2) {
//             bitvector<N> neighbours = g.out(u) | g.in(u);
//             visit(neighbours, [&](size_t v) {
//                 if (intersection_size(g.out(v) & g.in(v), neighbours) > 0)
//                     res = true;
//             });
//             // 2 out, double edge between them gives reduction
//         } else if (g.degree(u) == 2) {
//             size_t n[2];
//             visit(g.out(u), [&, i = 0](size_t v) mutable { n[i++] = v; });
//             if (test(g.out(n[0]), n[1]) && test(g.out(n[1]), n[0])) {
//                 res = true;
//             }
//             // 2 in, double edge between them gives reduction
//         } else if (g.in_degree(u) == 2) {
//             size_t n[2], i = 0;
//             visit(g.in(u), [&](size_t v) { n[i++] = v; });
//             if (test(g.out(n[0]), n[1]) && test(g.out(n[1]), n[0])) {
//                 res = true;
//             }
//         }
//         if (res) {
//             exclude_from_fvs(g, fvs, gs, u);
//             return true;
//         }
//     }

//     return false;
// }

// bool redundant_edge_meta(graph &g, bitvector<N> &fvs, graph_search &gs, size_t u) {
//     assert(u < N * 64 && test(g.active_vertices(), u));

//     std::function<bool(bitvector<N> &, size_t, bool)> test_walk = [&](bitvector<N> &path, size_t v, size_t d) {
//         while (true) {
//             if (intersection_size(path, g.out(v)) > 0)
//                 return true;
//             bitvector<N> n = g.out(v) & ~g.in(v);
//             if (test(n, u) || (d >= 1 && popcount(n) != 1))
//                 return false;
//             if (popcount(n) != 1) {
//                 bitvector<N> branched_path;
//                 bool res = true;
//                 visit(n, [&](size_t w) {
//                     branched_path = path;
//                     set(branched_path, w);
//                     res &= test_walk(branched_path, w, d + 1);
//                 });
//                 return res;
//             }
//             v = first(n);
//             set(path, v);
//         }
//     };

//     bool res = false;
//     bitvector<N> path{};
//     visit(g.out(u) & ~g.in(u), [&](size_t v) {
//         path.fill(0ull);
//         set(path, v);
//         if (test_walk(path, v, 0)) {
//             res = true;
//             g.remove_edge(u, v);
//             push_search(gs, v);
//         }
//     });
//     if (res) {
//         push_search(gs, u);
//     }

//     return res;
// }

bool SCC_edge_reduction(graph &g, bitvector &fvs, const bitvector &nodes, graph_search &gs, std::vector<bitvector> &SCC) {
    gs.DFS_stack.clear();
    gs.L.clear();
    gs.DFS_visited.clear();
    SCC.clear();

    std::function<void(size_t)> dfs_visit = [&](size_t u) {
        if (gs.DFS_visited.get(u))
            return;
        gs.DFS_visited.set(u);
        for (size_t v : g.out(u))
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
        for (size_t v : g.in(u))
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

// bool SCC_edge_reduction_L2(graph &g, bitvector<N> &fvs, const bitvector<N> &nodes, graph_search &gs) {
//     gs.DFS_stack.clear();
//     gs.L.clear();
//     gs.DFS_visited.fill(0ull);

//     std::function<void(size_t)> dfs_visit = [&](size_t u) {
//         if (test(gs.DFS_visited, u))
//             return;
//         set(gs.DFS_visited, u);
//         visit(g.out(u) & ~g.in(u), [&](size_t v) { dfs_visit(v); });
//         gs.L.push_back(u);
//     };

//     visit(nodes, [&](size_t u) {
//         gs.SCC_id[u] = N * 64;
//         dfs_visit(u);
//     });

//     std::function<void(size_t, size_t)> assign = [&](size_t u, size_t root) {
//         if (gs.SCC_id[u] != N * 64)
//             return;
//         gs.SCC_id[u] = root;
//         visit(g.in(u) & ~g.out(u), [&](size_t v) { assign(v, root); });
//     };

//     for (auto u : gs.L) {
//         if (gs.SCC_id[u] != N * 64)
//             continue;
//         assign(u, u);
//     }

//     bool res = false;
//     visit(nodes, [&](size_t u) {
//         visit(g.out(u) & ~g.in(u), [&](size_t v) {
//             if (gs.SCC_id[u] != gs.SCC_id[v]) {
//                 g.remove_edge(u, v);
//                 push_search(gs, u);
//                 push_search(gs, v);
//                 res = true;
//             }
//         });
//     });

//     return res;
// }

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

void reduce_graph(graph &g, bitvector &fvs, const bitvector &nodes, graph_search &gs, std::vector<bitvector> &SCC) {
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
            case reductions::twin_vertices:
                found = twin_vertices_reduction(g, fvs, gs, u);
                break;
            case reductions::clique_and_one:
                found = clique_and_one(g, fvs, gs, u);
                break;
            case reductions::square_fold:
                found = square_fold(g, fvs, gs, u);
                break;
            // case reductions::specific_pattern:
            //     found = specific_pattern(g, fvs, gs, u);
            //     break;
            // case reductions::redundant_edge_meta:
            //     found = redundant_edge_meta(g, fvs, gs, u);
            //     break;
            default:
                break;
            }
            if (found)
                rule = 0;
        }

        if (rule == num_reductions) {
            gs.tmp2.set_and(nodes, g.active_vertices());
            if (SCC_edge_reduction(g, fvs, gs.tmp2, gs, SCC))
                rule = 0;
        }
        // if (rule == num_reductions && SCC_edge_reduction_L2(g, fvs, nodes & g.active_vertices(), gs)) {
        //     rule = 0;
        // }
    }
}
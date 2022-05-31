#include "bounds.hpp"
#include "local_search.hpp"
#include "reduction_engine.hpp"
#include "reductions.hpp"
#include <algorithm>
#include <numeric>

bool has_cycle(const sparse_graph &g, size_t u, std::vector<uint32_t> &prev) {
    std::vector<uint32_t> current, next;
    bitvector visited(g.size());

    visited.clear();
    current.clear();
    current.push_back(u);
    next.clear();
    std::fill(std::begin(prev), std::end(prev), g.size());

    while (current.size() > 0) {
        next.clear();
        for (auto v : current) {
            for (auto w : g.out(v)) {
                if (visited.get(w))
                    continue;
                if (w == u) {
                    prev[u] = v;
                    return true;
                }
                prev[w] = v;
                next.push_back(w);
                visited.set(w);
            }
        }
        std::swap(next, current);
    }
    return false;
}

bool has_length_three_cycle(const sparse_graph &g, size_t u) {
    if (!g.is_active(u))
        return false;
    for (auto in : g.in(u)) {
        for (auto out : g.out(u)) {
            if (g.has_edge(out, in))
                return true;
        }
    }
    return false;
}

uint32_t reducing_peeling_lower_bound(sparse_graph &g) {
    bitvector fvs(g.size());
    reduction_engine re;
    graph_search gs(g.size(), false);

    std::vector<uint32_t> vertices(g.size()), prev(g.size());
    std::iota(std::begin(vertices), std::end(vertices), 0);
    std::sort(std::begin(vertices), std::end(vertices), [&](auto a, auto b) {
        if (g.pi_degree(a) > 0 || g.pi_degree(b) > 0)
            return g.pi_degree(a) > g.pi_degree(b);
        bool ac = has_length_three_cycle(g, a), bc = has_length_three_cycle(g, b);
        return (ac && !bc) || (ac == bc && g.out_degree(a) + g.in_degree(a) > g.out_degree(b) + g.in_degree(b));
    });

    for (auto u : vertices) {
        if (!g.is_active(u))
            continue;
        if (g.pi_degree(u) > 0) { // length-2 clique
            gs.tmp.clear();
            for (auto v : g.pi(u)) {
                if (std::all_of(std::begin(gs.tmp), std::end(gs.tmp), [&](auto w) { return g.has_edge(v, w) && g.has_edge(w, v); })) {
                    gs.tmp.push_back(v);
                }
            }
            for (auto v : gs.tmp) {
                add_to_fvs(g, fvs, gs, re, v);
            }
            queue_neighbourhood(g, gs, u);
            re.remove_include_vertex(g, u);
        } else { // find cycle
            if (has_cycle(g, u, prev)) {
                uint32_t v = prev[u];
                while (v != u) {
                    queue_neighbourhood(g, gs, v);
                    re.remove_include_vertex(g, v);
                    v = prev[v];
                }
                add_to_fvs(g, fvs, gs, re, u);
            } else {
                exclude_from_fvs(g, fvs, gs, re, u);
            }
        }
        reduce_graph(g, fvs, gs, re);
    }

    re.unfold_graph(g, 0, fvs);

    return fvs.popcount();
}

bitvector local_search_upper_bound(const sparse_graph &g, size_t iterations) {
    local_search ls(g, 0.25, 0);
    ls.greedy_one_zero_swaps();

    double T = 0.25;

    while (iterations > 0) {
        ls.set_temperature(T);
        ls.search(g, g.size() * (1.0 - ((T - 0.05) / (0.35 - 0.05))));

        T *= 0.999;
        if (T < 0.05) {
            ls.greedy_one_zero_swaps();
            ls.greedy_one_zero_swaps_dfs(g);
            if (ls.get_current_cost() > ls.get_best_cost())
                ls.return_to_best(g);

            T = 0.25;
            iterations--;
        }
    }

    return ls.get_best();
}
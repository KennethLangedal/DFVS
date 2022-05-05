#include "solver.hpp"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <optional>
#include <sstream>

bool has_cycle(const graph &g, const bitvector &tabu, size_t u, size_t &length, std::vector<size_t> &prev) {
    static bitvector current(g.size()), next(g.size()), visited(g.size()), tmp(g.size());
    if (current.size() != g.size()) {
        current = bitvector(g.size()), next = bitvector(g.size()), visited = bitvector(g.size()), tmp = bitvector(g.size());
    }
    current.set_and_not(g.out(u), g.in(u));
    current.set_and_not(current, tabu);
    visited = tabu;
    visited |= current;
    std::fill(std::begin(prev), std::end(prev), g.size());
    for (size_t v : current)
        prev[v] = u;
    length = 0;
    while (current.popcount() > 0) {
        ++length;
        next.clear();
        if (current.get(u))
            return true;
        for (size_t v : current) {
            tmp.set_and_not(g.out(v), g.in(v));
            tmp.set_and_not(tmp, visited);
            for (size_t w : tmp) {
                if (prev[w] == g.size() || g.in_degree(prev[w]) + g.out_degree(prev[w]) > g.in_degree(v) + g.out_degree(v)) {
                    prev[w] = v;
                }
                next.set(w);
                visited.set(w);
            }
        }
        std::swap(next, current);
    }
    return false;
}

// size_t lower_bound_VC(const graph &g, const bitvector &nodes, const std::vector<size_t> &order, bitvector &visited) {
//     static std::vector<size_t> set_id(g.size()), set_sizes(g.size()), set_counts(g.size());
//     static std::vector<bool> set_full(g.size());
//     if (set_id.size() != g.size()) {
//         set_id.resize(g.size());
//         set_sizes.resize(g.size());
//         set_counts.resize(g.size());
//         set_full.resize(g.size());
//     }

//     size_t res = 0;
//     std::fill(std::begin(set_sizes), std::end(set_sizes), 0);
//     for (size_t u : nodes) {
//         set_id[u] = g.size();
//     }
//     for (size_t u : order) {

//     }

//     return 0;
// }

// Lower bound counting length 2 cycle cliques and general cycles
size_t lower_bound(const graph &g, const bitvector &nodes) {
    static bitvector visited(g.size()), tmp(g.size()), tmp2(g.size()), tmp3(g.size()), tmp4(g.size());
    static std::vector<size_t> order;
    if (visited.size() != g.size()) {
        visited = bitvector(g.size()), tmp = bitvector(g.size()), tmp2 = bitvector(g.size()), tmp3 = bitvector(g.size()), tmp4 = bitvector(g.size());
    }
    order.clear();
    for (size_t u : nodes) {
        order.push_back(u);
    }
    std::sort(std::begin(order), std::end(order), [&](size_t a, size_t b) { return g.out_degree(a) + g.in_degree(a) < g.out_degree(b) + g.in_degree(b); });
    visited.set_and_not(nodes, g.active_vertices());

    size_t res = 0;
    // length 2 cycle cliques (2 and 3)
    for (size_t u : order) {
        if (visited.get(u))
            continue;
        tmp.set_and(g.out(u), g.in(u));
        tmp.set_and_not(tmp, visited);
        if (tmp.popcount() > 0) {
            size_t v = *tmp.begin();
            if (tmp.popcount() > 1) {
                tmp2.set_and(g.out(v), g.in(v));
                tmp2.set_and_not(tmp2, visited);
                if (tmp2.intersection_size(tmp) > 0) {
                    tmp3.set_and(tmp, tmp2);
                    size_t w = *tmp3.begin();
                    visited.set(u);
                    visited.set(v);
                    visited.set(w);
                    res += 2;
                } else {
                    visited.set(u);
                    visited.set(v);
                    res++;
                }
            } else {
                visited.set(u);
                visited.set(v);
                res++;
            }
        }
    }

    static std::vector<size_t> prev(g.size());
    if (prev.size() != g.size()) {
        prev.resize(g.size());
    }
    size_t length;

    // General case
    for (size_t u : order) {
        if (visited.get(u))
            continue;
        if (has_cycle(g, visited, u, length, prev)) {
            size_t v = prev[u];
            while (v != u) {
                visited.set(v);
                v = prev[v];
            }
            visited.set(u);
            res++;
        }
    }

    return res;
}

size_t find_dense_branch_node(const graph &g, const bitvector &nodes) {
    size_t best = g.size(), best_cnt;
    for (size_t u : nodes) {
        size_t u_cnt = g.in_degree(u) + g.out_degree(u); // g.in_degree(u) + g.out_degree(u);
        if (best == g.size() || u_cnt > best_cnt) {
            best = u;
            best_cnt = u_cnt;
        }
    }

    return best;
}

bool clique(const graph &g, const bitvector &n) {
    if (n.popcount() < 2)
        return true;
    for (size_t u : n) {
        if (g.out(u).intersection_size(n) != n.popcount() - 1)
            return false;
    }
    return true;
}

// When deciding to include u
void add_or_exclude_extra(graph &g, bitvector &fvs, graph_search &gs, size_t u) {
    if (g.out(u).intersection_size(g.in(u)) == 0) {
        if (g.out(u).popcount() == 2) {
            for (size_t v : g.out(u)) {
                exclude_from_fvs(g, fvs, gs, v);
            }
        }
        if (g.in(u).popcount() == 2) {
            for (size_t v : g.in(u)) {
                exclude_from_fvs(g, fvs, gs, v);
            }
        }
    } else if (g.out(u) == g.in(u)) { // Mirrors
        gs.tmp.set_and(g.out(u), g.in(u));
        gs.tmp1.clear();
        for (size_t v : gs.tmp) {
            gs.tmp2.set_and(g.out(v), g.in(v));
            gs.tmp1 |= gs.tmp2;
        }
        gs.tmp1.set_and_not(gs.tmp1, gs.tmp);
        for (size_t v : gs.tmp1) {
            gs.tmp2.set_and(g.out(v), g.in(v));
            gs.tmp.set_and_not(gs.tmp, gs.tmp2);
            if (gs.tmp.popcount() == 0 || clique(g, gs.tmp)) { // empty or clique
                add_to_fvs(g, fvs, gs, v);
            }
            gs.tmp.set_and(g.out(u), g.in(u));
        }
    } else if (g.out(u).intersection_size(g.in(u)) == 1) {
        if (g.out_degree(u) == 2) {
            gs.tmp.set_and_not(g.out(u), g.in(u));
            size_t v = *gs.tmp.begin();
            exclude_from_fvs(g, fvs, gs, v);
        }
        if (g.in_degree(u) == 2) {
            gs.tmp.set_and_not(g.in(u), g.out(u));
            size_t v = *gs.tmp.begin();
            exclude_from_fvs(g, fvs, gs, v);
        }
    }
}

bitvector solve_req(graph &g, graph_search &gs, size_t cost, size_t ub, size_t d, size_t &lb_counter) {
    size_t t0 = g.get_timestamp();
    bitvector fvs(g.size());

    std::vector<bitvector> SCC;
    reduce_graph(g, fvs, g.active_vertices(), gs, SCC);

    if (g.active_vertices().popcount() == 0) {
#ifdef VERBOSE
        std::cout << "\x1b[2K" << d << " " << ub << " " << lb_counter << '\r' << std::flush;
#endif
        g.unfold_graph(t0, fvs);
        return fvs;
    }

    if (cost + fvs.popcount() + lower_bound(g, g.active_vertices()) >= ub) {
        fvs.fill();
        lb_counter++;
#ifdef VERBOSE
        std::cout << "\x1b[2K" << d << " " << ub << " " << lb_counter << '\r' << std::flush;
#endif
        return fvs;
    }

    size_t u = find_dense_branch_node(g, g.active_vertices());
    size_t t1 = g.get_timestamp();

    // include u
    bitvector fvs_inc_u(g.size());
    add_or_exclude_extra(g, fvs_inc_u, gs, u);
    add_to_fvs(g, fvs_inc_u, gs, u);
    fvs_inc_u |= solve_req(g, gs, cost + fvs_inc_u.popcount() + fvs.popcount(), ub, d + 1, lb_counter);
    g.unfold_graph(t1, fvs_inc_u);
    ub = std::min(fvs_inc_u.popcount() + fvs.popcount() + cost, ub);

    // exclude u
    bitvector fvs_exc_u(g.size());
    exclude_from_fvs(g, fvs_exc_u, gs, u);
    fvs_exc_u |= solve_req(g, gs, cost + fvs_exc_u.popcount() + fvs.popcount(), ub, d + 1, lb_counter);
    g.unfold_graph(t1, fvs_exc_u);

    if (fvs_inc_u.popcount() < fvs_exc_u.popcount()) {
        fvs |= fvs_inc_u;
    } else {
        fvs |= fvs_exc_u;
    }

    g.unfold_graph(t0, fvs);

    return fvs;
}

size_t k_param(const graph &g) {
    std::vector<size_t> degrees(g.size(), 0);
    for (size_t u : g.active_vertices()) {
        degrees[u] = g.out(u).intersection_size(g.in(u));
    }
    std::sort(std::begin(degrees), std::end(degrees), std::greater<size_t>());
    size_t k = 0;
    while (k < g.size() && degrees[k] > k)
        k++;
    return k;
}

#include "local_search.hpp"

bitvector solve(graph &g) {
    graph_search gs(g.size());
    std::vector<bitvector> SCC;
    bitvector res(g.size());

    reduce_graph(g, res, g.active_vertices(), gs, SCC, true);

    if (g.active_vertices().popcount() == 0) {
        g.unfold_graph(0, res);
        return res;
    }

#ifdef VERBOSE
    std::cout << SCC.size() << std::endl;
#endif

    for (auto &&c : SCC) {

        std::vector<size_t> org_label(c.popcount()), new_label(g.size(), g.size());
        size_t i = 0;
        for (size_t u : c) {
            new_label[u] = i;
            org_label[i++] = u;
        }

        std::string new_graph = std::to_string(c.popcount()) + " 0 0\n";
        size_t edges = 0;
        for (size_t u : c) {
            for (size_t v : g.out(u)) {
                new_graph += std::to_string(new_label[v] + 1) + " ";
                if (!g.in(u).get(v))
                    edges++;
            }
            new_graph += "\n";
        }

        graph _g;
        std::stringstream ss(new_graph);
        ss >> _g;

#ifdef VERBOSE
        std::cout << _g.size() << " " << res.popcount() << std::endl;
#endif

#ifdef PRINT_COMPONENT
        std::ofstream fs("scripts/plot_data");
        _g.print_edgelist(fs, _g.active_vertices());
        exit(0);
#endif
        // bitvector _fvs(_g.size());
        // _fvs.fill();

        // for (size_t i = 0; i < 10; ++i) {
        //     local_search ls(_g.size(), i);
        //     ls.search(_g);
        //     if (ls.get_best().popcount() < _fvs.popcount()) {
        //         _fvs = ls.get_best();
        //     }
        // }

        local_search ls(_g.size(), i);
        ls.search(_g);

#ifdef VERBOSE
        std::cout << "\x1b[2K" << ls.get_best().popcount() + res.popcount() << std::endl;
#endif

        graph_search _gs(_g.size());
        size_t lb_counter = 0;

        bitvector _fvs = solve_req(_g, _gs, 0, ls.get_best().popcount(), 0, lb_counter);
        // bitvector _fvs = solve_req(_g, _gs, 0, _g.size() + 1, 0, lb_counter);

        if (ls.get_best().popcount() < _fvs.popcount()) {
            _fvs = ls.get_best();
        }

#ifdef VERBOSE
        std::cout << "\x1b[2K" << _g.size() << " " << ls.get_best().popcount() << " " << lower_bound(_g, _g.active_vertices()) << " " << _fvs.popcount() << std::endl;
#endif
        // for (size_t u : _fvs) {
        //     res.set(org_label[u]);
        // }

        for (size_t u : _fvs) {
            res.set(org_label[u]);
        }
    }

    g.unfold_graph(0, res);
    return res;
}
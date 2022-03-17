#include "solver.hpp"
#include <chrono>
#include <fstream>
#include <optional>

// Lower bound counting 2 and 3 length cycles
size_t lower_bound(const graph &g, const bitvector &nodes) {
    static bitvector visited(g.size()), tmp(g.size()), tmp2(g.size()), tmp3(g.size());
    visited.set_and_not(nodes, g.active_vertices());

    size_t res = 0;
    // length 2 cycles
    for (size_t u : nodes) {
        if (visited.get(u))
            continue;
        tmp.set_and(g.out(u), g.in(u));
        tmp2.set_and_not(tmp, visited);
        if (tmp2.popcount() > 0) {
            size_t v = *tmp2.begin();
            visited.set(u);
            visited.set(v);
            res++;
        }
    }

    // length 3 cycles
    for (size_t u : nodes) {
        if (visited.get(u))
            continue;
        tmp.set_and_not(g.out(u), visited);
        for (size_t v : tmp) {
            tmp2.set_and_not(g.out(v), visited);
            if (tmp2.intersection_size(tmp) > 0) {
                tmp3.set_and(tmp, tmp2);
                size_t w = *tmp3.begin();
                visited.set(u);
                visited.set(v);
                visited.set(w);
                res++;
                break;
            }
        }
    }

    return res;
}

size_t find_dense_branch_node(const graph &g, const bitvector &nodes) {
    size_t best = g.size(), best_cnt;
    for (size_t u : nodes) {
        size_t u_cnt = g.out_degree(u) + g.in_degree(u);
        if (u_cnt > 0 && (best == g.size() || u_cnt > best_cnt)) {
            best = u;
            best_cnt = u_cnt;
        }
    }

    return best;
}

size_t find_dense_branch_costly(graph &g, const bitvector &nodes, graph_search &gs) {
    size_t best = g.size(), best_v;
    bitvector fvs(g.size());
    fvs.set_and(g.active_vertices(), nodes);
    size_t n = fvs.popcount(), t = g.get_timestamp();
    fvs.clear();
    std::vector<bitvector> SCC;

    auto score = [&]() -> size_t {
        if (SCC.size() == 0)
            return 0;
        size_t res = g.size();
        for (auto cc : SCC) {
            if (cc.popcount() < res)
                res = cc.popcount();
        }
        return res;
    };

    for (size_t u : nodes) {
        if (!g.active_vertices().get(u))
            continue;

        add_to_fvs(g, fvs, gs, u);
        reduce_graph(g, fvs, nodes, gs, SCC);
        size_t score_inc = score();
        g.unfold_graph(t, fvs);
        fvs.clear();

        exclude_from_fvs(g, fvs, gs, u);
        reduce_graph(g, fvs, nodes, gs, SCC);
        size_t score_exc = score();
        g.unfold_graph(t, fvs);
        fvs.clear();

        if (std::max(score_inc, score_exc) < best) {
            best = std::max(score_inc, score_exc);
            best_v = u;
        }
    }

    return best_v;
}

static auto start_time = std::chrono::high_resolution_clock::now();

bitvector solve(graph &g, const bitvector &nodes, graph_search &gs, size_t cost, size_t ub, size_t d, size_t &lb_counter) {
    size_t t0 = g.get_timestamp();
    bitvector fvs(g.size());
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start_time);
    if (seconds.count() > 60) {
        std::cout << 60 << std::endl;
        exit(0);
    }
    std::vector<bitvector> SCC;
    reduce_graph(g, fvs, nodes, gs, SCC);
    if (nodes.intersection_size(g.active_vertices()) == 0) {
        g.unfold_graph(t0, fvs);
        return fvs;
    }
    // if (nodes.intersection_size(g.active_vertices()) < 30) {
    //     gs.tmp.set_and(nodes, g.active_vertices());
    //     bool test = true;
    //     // for (size_t u : gs.tmp) {
    //     //     if (g.out(u).intersection_size(g.in(u)) == 0)
    //     //         test = true;
    //     // }
    //     if (test) {
    //         std::cout << gs.tmp.popcount() << std::endl;
    //         std::ofstream fs("scripts/plot_data");
    //         g.print_edgelist(fs, gs.tmp);
    //         exit(0);
    //     }
    // }
    if (cost + fvs.popcount() + lower_bound(g, nodes) >= ub) {
        fvs.fill();
        lb_counter++;
        return fvs;
    }

    if (SCC.size() == 1) {
        size_t u = find_dense_branch_node(g, SCC.front());
        size_t t1 = g.get_timestamp();

        // exclude u
        bitvector fvs_ex_u(g.size());
        exclude_from_fvs(g, fvs_ex_u, gs, u);
        fvs_ex_u |= solve(g, SCC.front(), gs, cost + fvs.popcount() + fvs_ex_u.popcount(), ub, d + 1, lb_counter);
        g.unfold_graph(t1, fvs_ex_u);
        if (ub > cost + fvs.popcount() + fvs_ex_u.popcount())
            ub = cost + fvs.popcount() + fvs_ex_u.popcount();

        // include u
        bitvector fvs_inc_u(g.size());
        add_to_fvs(g, fvs_inc_u, gs, u);
        fvs_inc_u |= solve(g, SCC.front(), gs, cost + fvs.popcount() + fvs_inc_u.popcount(), ub, d + 1, lb_counter);
        g.unfold_graph(t1, fvs_inc_u);

        if (fvs_inc_u.popcount() < fvs_ex_u.popcount()) {
            fvs |= fvs_inc_u;
        } else {
            fvs |= fvs_ex_u;
        }
    } else {
        for (auto &&CC : SCC) {
            if (CC.popcount() == 0) {
                continue;
            }
            size_t u = find_dense_branch_node(g, CC);
            size_t t1 = g.get_timestamp();

            // exclude u
            bitvector fvs_ex_u(g.size());
            exclude_from_fvs(g, fvs_ex_u, gs, u);
            fvs_ex_u |= solve(g, CC, gs, fvs_ex_u.popcount(), g.size(), d + 1, lb_counter);
            g.unfold_graph(t1, fvs_ex_u);

            // include u
            bitvector fvs_inc_u(g.size());
            add_to_fvs(g, fvs_inc_u, gs, u);
            fvs_inc_u |= solve(g, CC, gs, fvs_inc_u.popcount(), fvs_ex_u.popcount(), d + 1, lb_counter);
            g.unfold_graph(t1, fvs_inc_u);

            if (fvs_inc_u.popcount() < fvs_ex_u.popcount()) {
                fvs |= fvs_inc_u;
            } else {
                fvs |= fvs_ex_u;
            }
        }
    }

    g.unfold_graph(t0, fvs);

    return fvs;
}
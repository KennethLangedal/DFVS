#include "solver.hpp"
#include "bounds.hpp"
#include "reductions.hpp"
#include <algorithm>
#include <fstream>

uint32_t find_dense_branch_node(const sparse_graph &g) {
    uint32_t best = g.size(), best_cnt;
    for (auto u : g.active_vertices()) {
        uint32_t u_cnt = (g.pi_degree(u) * 100) + g.in_degree_non_pi(u) + g.out_degree_non_pi(u);
        if (best == g.size() || u_cnt > best_cnt) {
            best = u;
            best_cnt = u_cnt;
        }
    }

    return best;
}

bool clique(const sparse_graph &g, const std::vector<uint32_t> &n) {
    if (n.size() < 2)
        return true;
    for (auto u : n) {
        if (intersection_size(g.pi(u), n) != n.size() - 1)
            return false;
    }
    return true;
}

// When deciding to include u
void add_or_exclude_extra(sparse_graph &g, bitvector &fvs, graph_search &gs, reduction_engine &re, uint32_t u) {
    if (g.pi_degree(u) == 0) {
        if (g.out_degree(u) == 2) {
            gs.tmp.clear();
            std::copy(std::begin(g.out(u)), std::end(g.out(u)), std::back_inserter(gs.tmp));
            for (auto v : gs.tmp) {
                exclude_from_fvs(g, fvs, gs, re, v);
            }
        }
        if (g.in_degree(u) == 2) {
            gs.tmp.clear();
            std::copy(std::begin(g.in(u)), std::end(g.in(u)), std::back_inserter(gs.tmp));
            for (auto v : gs.tmp) {
                exclude_from_fvs(g, fvs, gs, re, v);
            }
        }
    } else if (g.pi_only(u)) { // Mirrors
        std::vector<uint32_t> N2;
        for (auto v : g.pi(u)) {
            if (!g.pi_only(v))
                continue;
            gs.tmp.clear();
            std::set_union(std::begin(g.pi(v)), std::end(g.pi(v)), std::begin(N2), std::end(N2), std::back_inserter(gs.tmp));
            std::swap(N2, gs.tmp);
        }
        gs.tmp.clear();
        std::set_difference(std::begin(N2), std::end(N2), std::begin(g.pi(u)), std::end(g.pi(u)), std::back_inserter(gs.tmp));
        std::swap(N2, gs.tmp);
        for (auto v : N2) {
            if (!g.pi_only(v) || v == u)
                continue;
            gs.tmp.clear();
            std::set_difference(std::begin(g.pi(u)), std::end(g.pi(u)), std::begin(g.pi(v)), std::end(g.pi(v)), std::back_inserter(gs.tmp));
            if (gs.tmp.size() == 0 || clique(g, gs.tmp)) { // empty or clique
                add_to_fvs(g, fvs, gs, re, v);
            }
        }
    }
}

bitvector solve_req(sparse_graph &g, graph_search &gs, reduction_engine &re, uint32_t cost, size_t ub, uint32_t d, uint32_t &lb_counter) {
    size_t t0 = re.get_timestamp();
    bitvector fvs(g.size());

    reduce_graph(g, fvs, gs, re, true);

    if (g.active_vertices().popcount() == 0) {
#ifdef VERBOSE
        std::cout << "\x1b[2K" << d << " " << ub << " " << lb_counter << '\r' << std::flush;
#endif
        re.unfold_graph(g, t0, fvs);
        return fvs;
    }

    if (cost + fvs.popcount() + reducing_peeling_lower_bound(g) >= ub) {
        fvs.fill();
        lb_counter++;
#ifdef VERBOSE
        std::cout << "\x1b[2K" << d << " " << ub << " " << lb_counter << '\r' << std::flush;
#endif
        return fvs;
    }

    if (gs.SCC.size() > 1) {
        for (auto &&c : gs.SCC) {

            sparse_graph _g(g, c);
            graph_search _gs(_g.size(), false);

            bitvector _fvs = solve_req(_g, _gs, re, 0, _g.size(), d, lb_counter);

            for (auto u : _fvs) {
                fvs.set(_g.original_label(u));
            }
        }
    } else {
        uint32_t u = find_dense_branch_node(g);
        size_t t1 = re.get_timestamp();

        // include u
        bitvector fvs_inc_u(g.size());
        add_or_exclude_extra(g, fvs_inc_u, gs, re, u);
        add_to_fvs(g, fvs_inc_u, gs, re, u);
        fvs_inc_u |= solve_req(g, gs, re, cost + fvs_inc_u.popcount() + fvs.popcount(), ub, d + 1, lb_counter);
        re.unfold_graph(g, t1, fvs_inc_u);
        ub = std::min(fvs_inc_u.popcount() + fvs.popcount() + cost, ub);

        // exclude u
        bitvector fvs_exc_u(g.size());
        exclude_from_fvs(g, fvs_exc_u, gs, re, u);
        fvs_exc_u |= solve_req(g, gs, re, cost + fvs_exc_u.popcount() + fvs.popcount(), ub, d + 1, lb_counter);
        re.unfold_graph(g, t1, fvs_exc_u);

        if (fvs_inc_u.popcount() < fvs_exc_u.popcount()) {
            fvs |= fvs_inc_u;
        } else {
            fvs |= fvs_exc_u;
        }
    }

    re.unfold_graph(g, t0, fvs);

    return fvs;
}

bitvector solve(sparse_graph &g) {
    graph_search gs(g.size());
    reduction_engine re;
    bitvector fvs(g.size());

    reduce_graph(g, fvs, gs, re, true);

#ifdef VERBOSE
    std::cout << "Size before: " << g.size() << ", Size after: " << g.active_vertices().popcount() << ", Components: " << gs.SCC.size() << ", Res offset: " << fvs.popcount() << std::endl;
#endif

    if (g.active_vertices().popcount() == 0) {
        re.unfold_graph(g, 0, fvs);
        return fvs;
    }

    for (auto &&c : gs.SCC) {

        sparse_graph _g(g, c);

        bitvector ub = local_search_upper_bound(_g);

#ifdef VERBOSE
        std::cout << "Component size: " << _g.size() << ", LB " << reducing_peeling_lower_bound(_g) << ", UB: " << ub.popcount() << std::endl;

        size_t max_in_out = 0, non_pi = 0, pi = 0;
        for (auto u : _g.active_vertices()) {
            max_in_out = std::max(max_in_out, std::min(_g.out_degree(u), _g.in_degree(u)));
            non_pi += _g.out_degree_non_pi(u);
            pi += _g.pi_degree(u);
        }

        std::cout << "Heigest min degree " << max_in_out << ", Pi edges: " << pi << ", Non-pi edges: " << non_pi << std::endl;

        if (true || (non_pi > 0 && non_pi < 10 && _g.size() < 200)) {
            std::cout << "Printing graph" << std::endl;
            std::ofstream fs("scripts/plot_data");
            for (auto u : _g.active_vertices()) {
                for (auto v : _g.out(u)) {
                    fs << u << " " << v << std::endl;
                }
            }
            exit(0);
        }
#endif

        graph_search _gs(_g.size(), false);
        uint32_t lb_counter = 0;

        bitvector _fvs = solve_req(_g, _gs, re, 0, ub.popcount(), 0, lb_counter);

        if (ub.popcount() < _fvs.popcount()) {
            _fvs = ub;
        }

        for (auto u : _fvs) {
            fvs.set(_g.original_label(u));
        }
    }

    re.unfold_graph(g, 0, fvs);
    return fvs;
}
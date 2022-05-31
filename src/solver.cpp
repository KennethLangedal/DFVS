#include "solver.hpp"
#include "bounds.hpp"
#include "reductions.hpp"
#include <algorithm>
#include <fstream>

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

uint32_t find_dense_branch_node(const sparse_graph &g) {
    uint32_t best = g.size(), best_cnt;
    for (auto u : g.active_vertices()) {
        uint32_t u_cnt = g.in_degree(u) * g.out_degree(u); // (g.pi_degree(u) * 100) + (g.in_degree_non_pi(u) * g.out_degree_non_pi(u));
        if (best == g.size() || u_cnt > best_cnt) {
            best = u;
            best_cnt = u_cnt;
        }
    }

    return best;
}

std::vector<uint32_t> find_best_cycle_branch(sparse_graph &g, graph_search &gs, reduction_engine &re) {
    std::vector<uint32_t> res, size_after_include(g.size(), g.size()), tmp;
    size_t t0 = re.get_timestamp();
    bitvector fvs(g.size());
    for (auto u : g.active_vertices()) {
        add_to_fvs(g, fvs, gs, re, u);
        add_or_exclude_extra(g, fvs, gs, re, u);
        reduce_graph(g, fvs, gs, re, true);
        size_after_include[u] = g.active_vertices().popcount();
        if (gs.SCC.size() > 0) {
            size_after_include[u] = std::max_element(std::begin(gs.SCC), std::end(gs.SCC), [&](auto a, auto b) { return a.popcount() < b.popcount(); })->popcount();
        }
        re.unfold_graph(g, t0, fvs);
    }

    uint32_t best_branch = g.size();
    for (auto u : g.active_vertices()) {
        for (auto v : g.pi(u)) {
            uint32_t reduced_size = std::max(size_after_include[u], size_after_include[v]);
            if (reduced_size < best_branch || (reduced_size == best_branch && size_after_include[u] < size_after_include[res[0]])) {
                best_branch = reduced_size;
                res = {u, v};
            }
        }
        for (auto v : g.out_non_pi(u)) {
            tmp.clear();
            std::set_intersection(std::begin(g.in_non_pi(u)), std::end(g.in_non_pi(u)), std::begin(g.out_non_pi(v)), std::end(g.out_non_pi(v)), std::back_inserter(tmp));
            for (auto w : tmp) {
                uint32_t reduced_size = std::max(size_after_include[u], std::max(size_after_include[v], size_after_include[w]));
                if (reduced_size < best_branch || (reduced_size == best_branch && size_after_include[u] < size_after_include[res[0]])) {
                    best_branch = reduced_size;
                    res = {u, v, w};
                }
            }
        }
    }

#ifdef VERBOSE
    std::cout << "After best branch " << best_branch << " Lower Bound " << reducing_peeling_lower_bound(g) << std::endl;
#endif

    return res;
}

uint32_t find_best_branch_node(sparse_graph &g, graph_search &gs, reduction_engine &re) {
    bitvector fvs(g.size());

    uint32_t best_branch = g.size();
    uint32_t best_u;
    size_t t0 = re.get_timestamp();
    for (auto u : g.active_vertices()) {
        uint32_t inc, exc;

        add_to_fvs(g, fvs, gs, re, u);
        add_or_exclude_extra(g, fvs, gs, re, u);
        reduce_graph(g, fvs, gs, re, true);
        inc = g.active_vertices().popcount();
        if (gs.SCC.size() > 0) {
            inc = std::max_element(std::begin(gs.SCC), std::end(gs.SCC), [&](auto a, auto b) { return a.popcount() < b.popcount(); })->popcount();
        }
        re.unfold_graph(g, t0, fvs);

        exclude_from_fvs(g, fvs, gs, re, u);
        reduce_graph(g, fvs, gs, re, true);
        exc = g.active_vertices().popcount();
        if (gs.SCC.size() > 0) {
            exc = std::max_element(std::begin(gs.SCC), std::end(gs.SCC), [&](auto a, auto b) { return a.popcount() < b.popcount(); })->popcount();
        }
        re.unfold_graph(g, t0, fvs);

        if (std::max(inc, exc) < best_branch) {
            best_branch = std::max(inc, exc);
            best_u = u;
        }
    }

#ifdef VERBOSE
    std::cout << "After best branch " << best_branch << std::endl;
#endif

    return best_u;
}

bool update_constraints(sparse_graph &g, graph_search &gs, reduction_engine &re, bitvector fvs, std::vector<std::vector<uint32_t>> &c_v, std::vector<int32_t> &c_c) {
    for (size_t i = 0; i < c_v.size(); ++i) {
        for (auto u : c_v[i]) {
            if (fvs.get(u)) {
                c_c[i]--;
            }
        }
        if (c_c[i] == 0) {
            for (auto u : c_v[i]) {
                if (g.is_active(u)) {
                    exclude_from_fvs(g, fvs, gs, re, u);
                }
            }
        } else if (c_c[i] < 0) {
            return true;
        }
    }
    return false;
}

bitvector solve_req(sparse_graph &g, graph_search &gs, reduction_engine &re, uint32_t cost, size_t ub, uint32_t d, std::vector<std::vector<uint32_t>> &c_v, std::vector<int32_t> &c_c, uint32_t &lb_counter, std::string &log) {
    size_t t0 = re.get_timestamp();
    bitvector fvs(g.size());
    std::vector<int32_t> new_c_c = c_c;

    reduce_graph(g, fvs, gs, re, true);

    if (g.active_vertices().popcount() == 0) {
#ifdef VERBOSE
        std::cout << "\x1b[2K" << d << " " << ub << " " << lb_counter << " " << (log.size() < 50 ? log : "") << '\r' << std::flush;
#endif
        re.unfold_graph(g, t0, fvs);
        return fvs;
    }

    if (cost + fvs.popcount() + reducing_peeling_lower_bound(g) >= ub) {
        fvs.fill();
        lb_counter++;
#ifdef VERBOSE
        std::cout << "\x1b[2K" << d << " " << ub << " " << lb_counter << " " << (log.size() < 50 ? log : "") << '\r' << std::flush;
#endif
        return fvs;
    }

    if (gs.SCC.size() > 1) {
        for (auto &&c : gs.SCC) {

            sparse_graph _g(g, c);
            graph_search _gs(_g.size(), false);

            std::vector<std::vector<uint32_t>> _c_v;
            std::vector<int32_t> _c_c;

            bitvector _fvs = solve_req(_g, _gs, re, 0, _g.size(), d, _c_v, _c_c, lb_counter, log);

            for (auto u : _fvs) {
                fvs.set(_g.original_label(u));
            }
        }
    } else {
        uint32_t u = find_dense_branch_node(g);
        size_t t1 = re.get_timestamp();

        // include u
        log.push_back('I');
        bitvector fvs_inc_u(g.size());
        c_v.push_back(g.out(u));
        c_v.push_back(g.in(u));
        new_c_c.push_back(g.out_degree(u) - 2);
        new_c_c.push_back(g.in_degree(u) - 2);
        add_or_exclude_extra(g, fvs_inc_u, gs, re, u);
        add_to_fvs(g, fvs_inc_u, gs, re, u);
        if (!update_constraints(g, gs, re, fvs, c_v, new_c_c)) {
            fvs_inc_u |= solve_req(g, gs, re, cost + fvs_inc_u.popcount() + fvs.popcount(), ub, d + 1, c_v, new_c_c, lb_counter, log);
        } else {
            fvs_inc_u.fill();
        }
        c_v.pop_back();
        c_v.pop_back();
        re.unfold_graph(g, t1, fvs_inc_u);
        ub = std::min(fvs_inc_u.popcount() + fvs.popcount() + cost, ub);
        log.pop_back();

        // exclude u
        log.push_back('O');
        bitvector fvs_exc_u(g.size());
        exclude_from_fvs(g, fvs_exc_u, gs, re, u);
        fvs_exc_u |= solve_req(g, gs, re, cost + fvs_exc_u.popcount() + fvs.popcount(), ub, d + 1, c_v, new_c_c, lb_counter, log);
        re.unfold_graph(g, t1, fvs_exc_u);
        ub = std::min(fvs_exc_u.popcount() + fvs.popcount() + cost, ub);
        log.pop_back();

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
#endif

        graph_search _gs(_g.size(), false);
        uint32_t lb_counter = 0;

        std::string log;

        std::vector<std::vector<uint32_t>> c_v;
        std::vector<int32_t> c_c;

        bitvector _fvs = solve_req(_g, _gs, re, 0, ub.popcount(), 0, c_v, c_c, lb_counter, log);

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
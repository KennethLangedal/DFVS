#include "local_search.hpp"

size_t find_dense_node(const graph &g, const bitvector &nodes) {
    size_t best = g.size(), best_cnt;
    for (size_t u : nodes) {
        size_t u_cnt = g.in_degree(u) + g.out_degree(u);
        if (best == g.size() || u_cnt > best_cnt) {
            best = u;
            best_cnt = u_cnt;
        }
    }

    return best;
}

bool part_of_cycle(const graph &g, bitvector &visited, size_t u) {
    static bitvector current(g.size()), next(g.size()), tmp(g.size());
    current.set_and_not(g.out(u), visited);
    visited |= current;
    while (current.popcount() > 0) {
        next.clear();
        if (current.get(u))
            return true;
        for (size_t v : current) {
            next |= g.out(v);
        }
        next.set_and_not(next, visited);
        visited |= next;
        std::swap(next, current);
    }
    return false;
}

void reducing_peeling(graph &g, bitvector &dfvs, graph_search &gs, std::vector<bitvector> &SCC) {
    size_t t0 = g.get_timestamp();

    while (g.active_vertices().popcount() > 0) {
        size_t u = find_dense_node(g, g.active_vertices());
        add_to_fvs(g, dfvs, gs, u);
        reduce_graph(g, dfvs, g.active_vertices(), gs, SCC);
    }
    g.unfold_graph(t0, dfvs);
}

void remove_redundant(const graph &g, bitvector &dfvs) {
    bool found = true;
    static bitvector visited(g.size()), candidates(g.size());
    while (found) {
        found = false;
        candidates.set_and(g.active_vertices(), dfvs);
        for (size_t u : candidates) {
            visited = dfvs;
            visited.reset(u);
            if (!part_of_cycle(g, visited, u)) {
                found = true;
                dfvs.reset(u);
            }
        }
    }
}

void two_one_swap(const graph &g, bitvector &dfvs) {
    bool found = true;
    static bitvector candidates(g.size()), tmp(g.size()), tmp2(g.size());
    while (found) {
        found = false;
        candidates.set_and(g.active_vertices(), dfvs);
        for (size_t u : candidates) {
            if (!dfvs.get(u)) continue;

            size_t v = g.size();

            tmp.set_and_not(g.out(u), dfvs);
            

            if (tmp.popcount() == 1) {
                size_t v = *tmp.begin();
                tmp.set_or(g.out(v), g.in(v));
                tmp.set_and(tmp, dfvs);
                for (size_t w : tmp) {
                    if (w == u) continue;
                    tmp2.set_and_not(g.out(w), dfvs);
                    if (tmp2.popcount() == 1) {
                        found = true;
                        dfvs.reset(u);
                        dfvs.reset(w);
                        dfvs.set(v);
                        break;
                    }
                }
            }
            tmp.set_and_not(g.in(u), dfvs);
            if (tmp.popcount() == 1) {
                size_t v = *tmp.begin();
                tmp.set_and(g.out(v), dfvs);
                for (size_t w : tmp) {
                    if (w == u) continue;
                    tmp2.set_and_not(g.out(w), dfvs);
                    if (tmp2.popcount() == 1) {
                        found = true;
                        dfvs.reset(u);
                        dfvs.reset(w);
                        dfvs.set(v);
                        break;
                    }
                }
            }
        }
    }
}

#include <sstream>

void analysis(graph &g) {
    graph_search gs(g.size());
    std::vector<bitvector> SCC;
    bitvector res(g.size());
    reduce_graph(g, res, g.active_vertices(), gs, SCC);

    if (g.active_vertices().popcount() == 0) {
        g.unfold_graph(0, res);
        return;
    }

    std::vector<size_t> org_label(g.active_vertices().popcount()), new_label(g.size(), g.size());
    size_t i = 0;
    for (size_t u : g.active_vertices()) {
        new_label[u] = i;
        org_label[i++] = u;
    }

    std::string new_graph = std::to_string(g.active_vertices().popcount()) + " 0 0\n";
    size_t edges = 0;
    for (size_t u : g.active_vertices()) {
        for (size_t v : g.out(u)) {
            new_graph += std::to_string(new_label[v] + 1) + " ";
            if (!g.in(u).get(v))
                edges++;
        }
        new_graph += "\n";
    }

    std::cout << g.active_vertices().popcount() << " " << edges << std::endl;

    graph _g;
    std::stringstream ss(new_graph);
    ss >> _g;

    graph_search _gs(_g.size());
    size_t lb_counter = 0;
    bitvector _fvs(_g.size()), tmp;

    for (size_t u : _g.active_vertices()) {
        tmp = _fvs;
        add_to_fvs(_g, tmp, _gs, u);
        // reduce_graph(_g, tmp, _g.active_vertices(), _gs, SCC, false);
        _g.unfold_graph(0, tmp);
        tmp.reset(u);
        if (!part_of_cycle(_g, tmp, u)) {
            std::cout << "Found " << u << std::endl;
        }
    }
}
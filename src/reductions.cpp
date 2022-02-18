#include "reductions.hpp"

graph_search::graph_search() : search(num_reductions), visited(num_reductions, bitvector<N>{}) {
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

bool zero_degree_reduction(graph &g, std::vector<size_t> &fvs, graph_search &gs, size_t u) {
    assert(u < N * 64 && test(g.active_vertices(), u));
    if (g.in_degree(u) == 0 || g.degree(u) == 0) {
        deactivate_vertex(g, gs, u);
        return true;
    }
    return false;
}

bool self_edge_reduction(graph &g, std::vector<size_t> &fvs, graph_search &gs, size_t u) {
    assert(u < N * 64 && test(g.active_vertices(), u));
    if (g.self_loop(u)) {
        add_to_fvs(g, fvs, gs, u);
        return true;
    }
    return false;
}

bool one_degree_reduction(graph &g, std::vector<size_t> &fvs, graph_search &gs, size_t u) {
    assert(u < N * 64 && test(g.active_vertices(), u));
    if (g.in_degree(u) == 1) {
        visit(g.in(u), [&](size_t v) { push_search(gs, v); });
        visit(g.out(u), [&](size_t v) { push_search(gs, v); });
        g.deactive_single_in(u);
        return true;
    }
    if (g.degree(u) == 1) {
        visit(g.out(u), [&](size_t v) { push_search(gs, v); });
        visit(g.in(u), [&](size_t v) { push_search(gs, v); });
        g.deactive_single_out(u);
        return true;
    }
    return false;
}

bool isolated_vertex_reduction(graph &g, std::vector<size_t> &fvs, graph_search &gs, size_t u) {
    assert(u < N * 64 && test(g.active_vertices(), u));
    bitvector<N> tmp = g.out(u) & g.in(u);
    if (tmp == g.out(u) || tmp == g.in(u)) {
        bool res = true;
        size_t degree = popcount(tmp);
        visit(tmp, [&](size_t v) {
            res &= intersection_size(tmp, g.out(v)) == degree - 1;
        });
        if (res) {
            visit(tmp, [&](size_t v) { add_to_fvs(g, fvs, gs, v); });
            deactivate_vertex(g, gs, u);
            return true;
        }
    }
    return false;
}

bool twin_vertices_reduction(graph &g, std::vector<size_t> &fvs, graph_search &gs, size_t u) {
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

bool dominating_vertex_reduction(graph &g, std::vector<size_t> &fvs, graph_search &gs, size_t u) {
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

void add_to_fvs(graph &g, std::vector<size_t> &fvs, graph_search &gs, size_t u) {
    fvs.push_back(u);
    deactivate_vertex(g, gs, u);
}

void deactivate_vertex(graph &g, graph_search &gs, size_t u) {
    visit(g.out(u), [&](size_t v) { push_search(gs, v); });
    visit(g.in(u), [&](size_t v) { push_search(gs, v); });
    g.deactive_vertex(u);
}

void reduce_graph(graph &g, std::vector<size_t> &fvs, graph_search &gs) {
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
            case reductions::isolated_vertex:
                found = isolated_vertex_reduction(g, fvs, gs, u);
                break;
            case reductions::twin_vertices:
                found = twin_vertices_reduction(g, fvs, gs, u);
                break;
            case reductions::dominating_vertex:
                found = dominating_vertex_reduction(g, fvs, gs, u);
                break;
            default:
                break;
            }
            if (found)
                rule = 0;
        }
    }
}
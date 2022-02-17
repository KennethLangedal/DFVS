#include "reductions.hpp"

graph_search::graph_search() : search(num_reductions), visited(num_reductions, bitvector<N>{}) {
    for (size_t r = 0; r < num_reductions; ++r) {
        for (size_t u = 0; u < N * 64; ++u) {
            search[r].push_back(u);
        }
    }
}

void push_search(graph_search& gs, size_t u) {
    assert(u < N * 64);
    for (size_t r = 0; r < num_reductions; ++r) {
        if (test(gs.visited[r], u)) gs.search[r].push_back(u);
        reset(gs.visited[r], u);
    }
}

size_t pop_search(graph_search& gs, size_t rule) {
    assert(rule < num_reductions && !gs.search[rule].empty());
    size_t u = gs.search[rule].back();
    gs.search[rule].pop_back();
    set(gs.visited[rule], u);
    return u;
}

bool zero_degree_reduction(graph& g, std::vector<size_t>& fvs, graph_search& gs, size_t u) {
    assert(u < N * 64 && test(g.active_vertices(), u));
    if (g.in_degree(u) == 0 || g.degree(u) == 0) {
        deactivate_vertex(g, gs, u);
        return true;
    }
    return false;
}

bool self_edge_reduction(graph& g, std::vector<size_t>& fvs, graph_search& gs, size_t u) {
    assert(u < N * 64 && test(g.active_vertices(), u));
    if (g.self_loop(u)) {
        add_to_fvs(g, fvs, gs, u);
        return true;
    }
    return false;
}

bool one_degree_reduction(graph& g, std::vector<size_t>& fvs, graph_search& gs, size_t u) {
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

void add_to_fvs(graph& g, std::vector<size_t>& fvs, graph_search& gs, size_t u) {
    fvs.push_back(u);
    deactivate_vertex(g, gs, u);
}

void deactivate_vertex(graph& g, graph_search& gs, size_t u) {
    visit(g.out(u), [&](size_t v) { push_search(gs, v); });
    visit(g.in(u), [&](size_t v) { push_search(gs, v); });
    g.deactive_vertex(u);
}

void reduce_graph(graph& g, std::vector<size_t>& fvs, graph_search& gs) {
    size_t rule = 0;
    while (rule < num_reductions) {
        if (gs.search[rule].empty()) {
            rule++;
        } else {
            size_t u = pop_search(gs, rule);
            if (!test(g.active_vertices(), u)) continue;
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
                default:
                    break;
            }
            if (found) rule = 0;
        }
    }
}
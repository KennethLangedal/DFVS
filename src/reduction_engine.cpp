#include "reduction_engine.hpp"
#include <cassert>

void reduction_engine::add_edge(sparse_graph &g, uint32_t u, uint32_t v) {
    if (g.has_edge(u, v))
        return;
    g.add_edge(u, v);
    _log.push_back({reduction::add_edge, u, v, 0, 0, 0});
}

void reduction_engine::_undo_add_edge(sparse_graph &g, size_t u, size_t v) {
    assert(g.has_edge(u, v));
    g.remove_edge(u, v);
}

void reduction_engine::remove_edge(sparse_graph &g, uint32_t u, uint32_t v) {
    if (!g.has_edge(u, v))
        return;
    g.remove_edge(u, v);
    _log.push_back({reduction::remove_edge, u, v, 0, 0, 0});
}

void reduction_engine::_undo_remove_edge(sparse_graph &g, size_t u, size_t v) {
    assert(!g.has_edge(u, v));
    g.add_edge(u, v);
}

void reduction_engine::remove_include_vertex(sparse_graph &g, uint32_t u) {
    g.remove_vertex(u);
    _log.push_back({reduction::remove_vertex, u, 0, 0, 0, 0});
}

void reduction_engine::remove_exclude_vertex(sparse_graph &g, uint32_t u) {
    for (auto v : g.in_non_pi(u)) {
        for (auto w : g.out_non_pi(u)) {
            add_edge(g, v, w);
        }
    }
    for (auto v : g.pi(u)) {
        add_edge(g, v, v);
    }
    g.remove_vertex(u);
    _log.push_back({reduction::remove_vertex, u, 0, 0, 0, 0});
}

void reduction_engine::_undo_remove_vertex(sparse_graph &g, size_t u) {
    g.reactivate_vertex(u);
}

void reduction_engine::fold_twin(sparse_graph &g, uint32_t u, uint32_t v, bitvector &fvs) {
    uint32_t n[3], i = 0;
    for (auto x : g.pi(u))
        n[i++] = x;
    remove_include_vertex(g, u);
    remove_include_vertex(g, v);
    for (auto w : g.pi(n[1])) {
        add_edge(g, n[0], w);
        add_edge(g, w, n[0]);
    }
    for (auto w : g.pi(n[2])) {
        add_edge(g, n[0], w);
        add_edge(g, w, n[0]);
    }
    remove_include_vertex(g, n[1]);
    remove_include_vertex(g, n[2]);

    fvs.set(u);
    fvs.set(v);
    _log.push_back({reduction::fold_twin, u, v, n[0], 0, 0});
}

void reduction_engine::fold_clique_and_one(sparse_graph &g, uint32_t u, uint32_t v, bitvector &fvs) {
    for (auto w1 : g.out(u)) { // for every vertex in clique
        if (w1 == v)
            continue;
        for (auto w2 : g.out(v)) { // add v's edges
            if (w2 == w1)
                continue;
            add_edge(g, w1, w2);
        }
        for (auto w2 : g.in(v)) {
            if (w2 == w1)
                continue;
            add_edge(g, w2, w1);
        }
    }

    fvs.set(u);

    remove_include_vertex(g, u);
    remove_include_vertex(g, v);
    _log.push_back({reduction::fold_clique_and_one, u, v, 0, 0, 0});
}

void reduction_engine::fold_square(sparse_graph &g, uint32_t u, uint32_t v1, uint32_t v2, uint32_t w1, uint32_t w2, bitvector &fvs) {
    for (auto x : g.out(v2)) {
        if (x != v1)
            add_edge(g, v1, x);
    }
    for (auto x : g.in(v2)) {
        if (x != v1)
            add_edge(g, x, v1);
    }
    for (auto x : g.out(w2)) {
        if (x != w1)
            add_edge(g, w1, x);
    }
    for (auto x : g.in(w2)) {
        if (x != w1)
            add_edge(g, x, w1);
    }

    fvs.set(v2);
    fvs.set(w2);

    remove_include_vertex(g, u);
    remove_include_vertex(g, v2);
    remove_include_vertex(g, w2);
    _log.push_back({reduction::fold_square, u, v1, v2, w1, w2});
}

size_t reduction_engine::get_timestamp() const {
    return _log.size();
}

void reduction_engine::unfold_graph(sparse_graph &g, size_t time, bitvector &fvs) {
    while (_log.size() > time) {
        auto [t, u, v1, v2, w1, w2] = _log.back();
        bool all;
        _log.pop_back();
        switch (t) {
        case reduction::add_edge:
            _undo_add_edge(g, u, v1);
            break;
        case reduction::remove_edge:
            _undo_remove_edge(g, u, v1);
            break;
        case reduction::remove_vertex:
            _undo_remove_vertex(g, u);
            break;
        case reduction::fold_twin:
            if (fvs.get(v2)) {
                fvs.reset(u);
                fvs.reset(v1);

                for (auto x : g.pi(u)) {
                    fvs.set(x);
                }
            }
            break;
        case reduction::fold_clique_and_one:
            all = true;
            for (auto x : g.out(u)) {
                if (x != v1 && !fvs.get(x))
                    all = false;
            }
            if (all) {
                fvs.reset(u);
                fvs.set(v1);
            }
            break;
        case reduction::fold_square:
            fvs.reset(v2);
            fvs.reset(w2);
            if (fvs.get(v1)) {
                fvs.set(v2);
            }
            if (fvs.get(w1)) {
                fvs.set(w2);
            }
            if (!(fvs.get(v1) && fvs.get(w1))) {
                fvs.set(u);
            }
            break;

        default:
            break;
        }
    }
}
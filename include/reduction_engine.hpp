#pragma once
#include "bitvector.hpp"
#include <cassert>

enum class reduction {
    add_edge,
    remove_edge,
    remove_vertex,
    fold_clique_and_one,
    fold_square
};

struct action {
    reduction t;
    uint32_t u, v1, v2, w1, w2;
};

template <class graph>
class reduction_engine {
private:
    std::vector<action> _log;

    void _undo_add_edge(graph &g, size_t u, size_t v);
    void _undo_remove_edge(graph &g, size_t u, size_t v);
    void _undo_remove_vertex(graph &g, size_t u);

public:
    void add_edge(graph &g, uint32_t u, uint32_t v);
    void remove_edge(graph &g, uint32_t u, uint32_t v);
    void remove_include_vertex(graph &g, uint32_t u);
    void remove_exclude_vertex(graph &g, uint32_t u);
    void fold_clique_and_one(graph &g, uint32_t u, uint32_t v, bitvector &fvs);
    void fold_square(graph &g, uint32_t u, uint32_t v1, uint32_t v2, uint32_t w1, uint32_t w2, bitvector &fvs);

    size_t get_timestamp() const;

    void unfold_graph(graph &g, size_t time, bitvector &fvs);
};

template <class graph>
void reduction_engine<graph>::add_edge(graph &g, uint32_t u, uint32_t v) {
    if (g.has_edge(u, v))
        return;
    g.add_edge(u, v);
    _log.push_back({reduction::add_edge, u, v, 0, 0, 0});
}

template <class graph>
void reduction_engine<graph>::_undo_add_edge(graph &g, size_t u, size_t v) {
    assert(g.has_edge(u, v));
    g.remove_edge(u, v);
}

template <class graph>
void reduction_engine<graph>::remove_edge(graph &g, uint32_t u, uint32_t v) {
    if (!g.has_edge(u, v))
        return;
    g.remove_edge(u, v);
    _log.push_back({reduction::remove_edge, u, v, 0, 0, 0});
}

template <class graph>
void reduction_engine<graph>::_undo_remove_edge(graph &g, size_t u, size_t v) {
    assert(!g.has_edge(u, v));
    g.add_edge(u, v);
}

template <class graph>
void reduction_engine<graph>::remove_include_vertex(graph &g, uint32_t u) {
    g.remove_vertex(u);
    _log.push_back({reduction::remove_vertex, u, 0, 0, 0, 0});
}

template <class graph>
void reduction_engine<graph>::remove_exclude_vertex(graph &g, uint32_t u) {
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

template <class graph>
void reduction_engine<graph>::_undo_remove_vertex(graph &g, size_t u) {
    g.reactivate_vertex(u);
}

template <class graph>
void reduction_engine<graph>::fold_clique_and_one(graph &g, uint32_t u, uint32_t v, bitvector &fvs) {
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

template <class graph>
void reduction_engine<graph>::fold_square(graph &g, uint32_t u, uint32_t v1, uint32_t v2, uint32_t w1, uint32_t w2, bitvector &fvs) {
    for (auto x : g.out(v2)) {
        if (x != v1)
            add_edge(v1, x);
    }
    for (auto x : g.in(v2)) {
        if (x != v1)
            add_edge(x, v1);
    }
    for (auto x : g.out(w2)) {
        if (x != w1)
            add_edge(w1, x);
    }
    for (auto x : g.in(w2)) {
        if (x != w1)
            add_edge(x, w1);
    }

    fvs.set(v2);
    fvs.set(w2);

    remove_include_vertex(g, u);
    remove_include_vertex(g, v2);
    remove_include_vertex(g, w2);
    _log.push_back({reduction::fold_square, u, v1, v2, w1, w2});
}

template <class graph>
size_t reduction_engine<graph>::get_timestamp() const {
    return _log.size();
}

template <class graph>
void reduction_engine<graph>::unfold_graph(graph &g, size_t time, bitvector &fvs) {
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
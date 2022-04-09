#include "graph.hpp"
#include <cassert>

void graph::add_edge_raw(size_t u, size_t v) {
    _out_edges[u].set(v);
    _in_edges[v].set(u);
}

void graph::remove_edge_raw(size_t u, size_t v) {
    _out_edges[u].reset(v);
    _in_edges[v].reset(u);
}

void graph::add_edge(size_t u, size_t v) {
    assert(u < _N && v < _N && _active.get(u) && _active.get(v));
    if (_out_edges[u].get(v) || _in_edges[v].get(u))
        return;
    _log.push_back({action::add_edge, u, v, 0, 0, 0});
    _out_edges[u].set(v);
    _in_edges[v].set(u);
}

void graph::_undo_add_edge(size_t u, size_t v) {
    assert(u < _N && v < _N && _active.get(u) && _active.get(v) && _out_edges[u].get(v) && _in_edges[v].get(u));
    _out_edges[u].reset(v);
    _in_edges[v].reset(u);
}

void graph::remove_edge(size_t u, size_t v) {
    assert(u < _N && v < _N && _active.get(u) && _active.get(v) && _out_edges[u].get(v) && _in_edges[v].get(u));
    _log.push_back({action::remove_edge, u, v, 0, 0, 0});
    _out_edges[u].reset(v);
    _in_edges[v].reset(u);
}

void graph::_undo_remove_edge(size_t u, size_t v) {
    assert(u < _N && v < _N && _active.get(u) && _active.get(v) && !_out_edges[u].get(v) && !_in_edges[v].get(u));
    _out_edges[u].set(v);
    _in_edges[v].set(u);
}

void graph::remove_vertex(size_t u) {
    assert(u < _N && _active.get(u));
    _log.push_back({action::remove_vertex, u, 0, 0, 0, 0});
    _active.reset(u);
    for (size_t v : _out_edges[u])
        _in_edges[v].reset(u);
    for (size_t v : _in_edges[u])
        _out_edges[v].reset(u);
}

void graph::_undo_remove_vertex(size_t u) {
    assert(u < _N && !_active.get(u));
    _active.set(u);
    for (size_t v : _out_edges[u])
        _in_edges[v].set(u);
    for (size_t v : _in_edges[u])
        _out_edges[v].set(u);
}

void graph::fold_clique_and_one(size_t u, size_t v) {
    assert(u < _N && _active.get(u));

    for (size_t w1 : _out_edges[u]) {
        if (w1 == v)
            continue;
        for (size_t w2 : _out_edges[v]) {
            if (w2 == w1)
                continue;
            add_edge(w1, w2);
        }
        for (size_t w2 : _in_edges[v]) {
            if (w2 == w1)
                continue;
            add_edge(w2, w1);
        }
    }

    remove_vertex(u);
    remove_vertex(v);
    _log.push_back({action::fold_clique_and_one, u, v, 0, 0, 0});
}

void graph::fold_square(size_t u, size_t v1, size_t v2, size_t w1, size_t w2) {
    assert(u < _N && _out_edges[u].popcount() == 4 && _out_edges[u] == _in_edges[u]);
    for (size_t x : _out_edges[v2]) {
        if (x != v1)
            add_edge(v1, x);
    }
    for (size_t x : _in_edges[v2]) {
        if (x != v1)
            add_edge(x, v1);
    }
    for (size_t x : _out_edges[w2]) {
        if (x != w1)
            add_edge(w1, x);
    }
    for (size_t x : _in_edges[w2]) {
        if (x != w1)
            add_edge(x, w1);
    }
    remove_vertex(u);
    remove_vertex(v2);
    remove_vertex(w2);
    _log.push_back({action::fold_square, u, v1, v2, w1, w2});
}

size_t graph::out_degree(size_t u) const {
    assert(u < _N);
    return _out_edges[u].popcount();
}

size_t graph::in_degree(size_t u) const {
    assert(u < _N);
    return _in_edges[u].popcount();
}

bool graph::self_loop(size_t u) const {
    assert(u < _N);
    return _out_edges[u].get(u);
}

const bitvector &graph::active_vertices() const {
    return _active;
}

const bitvector &graph::out(size_t u) const {
    assert(u < _N);
    return _out_edges[u];
}

const bitvector &graph::in(size_t u) const {
    assert(u < _N);
    return _in_edges[u];
}

size_t graph::get_timestamp() const {
    return _log.size();
}

size_t graph::size() const {
    return _N;
}

void graph::unfold_graph(size_t time, bitvector &fvs) {
    while (_log.size() > time) {
        auto [t, u, v1, v2, w1, w2] = _log.back();
        bool all;
        _log.pop_back();
        switch (t) {
        case action::add_edge:
            _undo_add_edge(u, v1);
            break;
        case action::remove_edge:
            _undo_remove_edge(u, v1);
            break;
        case action::remove_vertex:
            _undo_remove_vertex(u);
            break;
        case action::fold_clique_and_one:
            all = true;
            for (size_t x : _out_edges[u]) {
                if (x != v1 && !fvs.get(x))
                    all = false;
            }
            if (all) {
                fvs.reset(u);
                fvs.set(v1);
            }
            break;
        case action::fold_square:
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

void graph::print_edgelist(std::ostream &os, const bitvector vertices) const {
    for (size_t u : vertices) {
        for (size_t v : _out_edges[u]) {
            os << u << " " << v << std::endl;
        }
    }
}

std::istream &operator>>(std::istream &is, graph &g) {
    size_t N_BUFFER = 1048576;
    std::vector<char> buffer(N_BUFFER);
    size_t buffer_i = N_BUFFER;

    auto get = [&](char &c) {
        if (buffer_i == N_BUFFER) {
            is.read(buffer.data(), N_BUFFER);
            buffer_i = 0;
        }
        c = buffer[buffer_i++];
    };

    auto peek = [&](char &c) {
        if (buffer_i == N_BUFFER) {
            is.read(buffer.data(), N_BUFFER);
            buffer_i = 0;
        }
        c = buffer[buffer_i];
    };

    char c;
    auto fscan = [&](auto &v) {
        get(c);
        v = 0;
        while (c >= '0' && c <= '9') {
            v = (v * 10) + (c - '0');
            get(c);
        }
    };

    peek(c);
    while (c == '%') {
        while (c != '\n')
            get(c);
        peek(c);
    }

    size_t N, E, tmp, u, v;
    fscan(N);
    fscan(E);
    fscan(tmp);

    g._N = N;
    g._active = bitvector(N);
    g._active.fill();
    g._out_edges.resize(N, bitvector(N));
    g._in_edges.resize(N, bitvector(N));

    u = 0;
    while (u < N) {
        fscan(v);
        if (v > 0) {
            g._out_edges[u].set(v - 1);
            g._in_edges[v - 1].set(u);
        }
        if (c == '\n') {
            u++;
            peek(c);
            while (c == '%') {
                while (c != '\n')
                    get(c);
                peek(c);
            }
        }
    }
    return is;
}
#include "dense_graph.hpp"
#include <algorithm>
#include <cassert>
#include <numeric>

size_t dense_graph::size() const {
    return _N;
}

bool dense_graph::is_active(uint32_t u) const {
    assert(u < _N);
    return _active.get(u);
}

bool dense_graph::has_edge(uint32_t u, uint32_t v) const {
    assert(is_active(u) && is_active(v));
    return _out_edges[u].get(v);
}

uint32_t dense_graph::original_label(uint32_t u) const {
    assert(u < _N);
    return _original_labels[u];
}

const bitvector &dense_graph::active_vertices() const {
    return _active;
}

void dense_graph::reactivate_vertex(uint32_t u) {
    assert(!is_active(u));

    for (auto v : in_non_pi(u)) {
        _out_edges_non_pi[v].set(u);
        _out_edges[v].set(u);
    }

    for (auto v : out_non_pi(u)) {
        _in_edges_non_pi[v].set(u);
        _in_edges[v].set(u);
    }

    for (auto v : pi(u)) {
        if (v == u)
            continue;
        _pi_edges[v].set(u);
        _out_edges[v].set(u);
        _in_edges[v].set(u);
    }

    _active.set(u);
}

void dense_graph::remove_vertex(uint32_t u) {
    assert(is_active(u));
    for (auto v : in_non_pi(u)) {
        _out_edges_non_pi[v].reset(u);
        _out_edges[v].reset(u);
    }

    for (auto v : out_non_pi(u)) {
        _in_edges_non_pi[v].reset(u);
        _in_edges[v].reset(u);
    }

    for (auto v : pi(u)) {
        if (v == u)
            continue;
        _pi_edges[v].reset(u);
        _out_edges[v].reset(u);
        _in_edges[v].reset(u);
    }

    _active.reset(u);
}

void dense_graph::remove_edge(uint32_t u, uint32_t v) {
    assert(is_active(u) && is_active(v) && has_edge(u, v));
    if (u == v) {
        _pi_edges[u].reset(u);
    } else if (has_edge(v, u)) {
        _pi_edges[u].reset(v);
        _pi_edges[v].reset(u);

        _out_edges_non_pi[v].set(u);
        _in_edges_non_pi[u].set(v);
    } else {
        _out_edges_non_pi[u].reset(v);
        _in_edges_non_pi[v].reset(u);
    }
    _out_edges[u].reset(v);
    _in_edges[v].reset(u);
}

void dense_graph::add_edge(uint32_t u, uint32_t v) {
    assert(u < _N && v < _N);
    if (has_edge(u, v))
        return;

    if (has_edge(v, u)) {
        _in_edges_non_pi[u].reset(v);
        _out_edges_non_pi[v].reset(u);

        _pi_edges[u].set(v);
        _pi_edges[v].set(u);
    } else if (u == v) {
        _pi_edges[u].set(u);
    } else {
        _out_edges_non_pi[u].set(v);
        _in_edges_non_pi[v].set(u);
    }
    _out_edges[u].set(v);
    _in_edges[v].set(u);
}

size_t dense_graph::out_degree(uint32_t u) const {
    assert(u < _N);
    return _out_edges[u].popcount();
}

size_t dense_graph::out_degree_non_pi(uint32_t u) const {
    assert(u < _N);
    return _out_edges_non_pi[u].popcount();
}

size_t dense_graph::in_degree(uint32_t u) const {
    assert(u < _N);
    return _in_edges[u].popcount();
}

size_t dense_graph::in_degree_non_pi(uint32_t u) const {
    assert(u < _N);
    return _in_edges_non_pi[u].popcount();
}

size_t dense_graph::pi_degree(uint32_t u) const {
    assert(u < _N);
    return _pi_edges[u].popcount();
}

const bitvector &dense_graph::out(uint32_t u) const {
    assert(u < _N);
    return _out_edges[u];
}

const bitvector &dense_graph::out_non_pi(uint32_t u) const {
    assert(u < _N);
    return _out_edges_non_pi[u];
}

const bitvector &dense_graph::in(uint32_t u) const {
    assert(u < _N);
    return _in_edges[u];
}

const bitvector &dense_graph::in_non_pi(uint32_t u) const {
    assert(u < _N);
    return _in_edges_non_pi[u];
}

const bitvector &dense_graph::pi(uint32_t u) const {
    assert(u < _N);
    return _pi_edges[u];
}

void dense_graph::parse_graph(std::istream &is, size_t N) {
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
    get(c);
    
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

    size_t u, v;

    _N = N;
    _active = bitvector(N);
    _active.fill();
    _out_edges.resize(N, bitvector(N));
    _out_edges_non_pi.resize(N, bitvector(N));
    _in_edges.resize(N, bitvector(N));
    _in_edges_non_pi.resize(N, bitvector(N));
    _pi_edges.resize(N, bitvector(N));
    _original_labels.resize(N);

    std::iota(std::begin(_original_labels), std::end(_original_labels), 0);

    u = 0;
    while (u < N) {
        fscan(v);
        if (v > 0) {
            _out_edges[u].set(v - 1);
            _in_edges[v - 1].set(u);
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

    for (uint32_t i = 0; i < N; ++i) {
        _pi_edges[i].set_and(_in_edges[i], _out_edges[i]);
        _in_edges_non_pi[i].set_and_not(_in_edges[i], _out_edges[i]);
        _out_edges_non_pi[i].set_and_not(_out_edges[i], _in_edges[i]);
    }
}

bool includes(const bitvector &a, const bitvector &b) {
    return b.subset_eq(a);
}

uint32_t intersection_size(const bitvector &a, const bitvector &b) {
    return a.intersection_size(b);
}
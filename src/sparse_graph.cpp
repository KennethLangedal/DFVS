#include "sparse_graph.hpp"
#include <algorithm>
#include <cassert>
#include <numeric>

sparse_graph::sparse_graph(const sparse_graph &g, const bitvector &nodes)
    : _N(nodes.popcount()), _active(_N),
      _out_edges(_N), _out_edges_non_pi(_N), _in_edges(_N), _in_edges_non_pi(_N), _pi_edges(_N),
      _original_labels(_N) {

    _active.fill();

    std::vector<uint32_t> new_labels(g.size());
    uint32_t i = 0;
    for (auto u : nodes) {
        new_labels[u] = i;
        _original_labels[i++] = u;
    }

    for (auto u : nodes) {
        for (auto v : g.out(u))
            _out_edges[new_labels[u]].push_back(new_labels[v]);
        for (auto v : g.out_non_pi(u))
            _out_edges_non_pi[new_labels[u]].push_back(new_labels[v]);
        for (auto v : g.in(u))
            _in_edges[new_labels[u]].push_back(new_labels[v]);
        for (auto v : g.in_non_pi(u))
            _in_edges_non_pi[new_labels[u]].push_back(new_labels[v]);
        for (auto v : g.pi(u))
            _pi_edges[new_labels[u]].push_back(new_labels[v]);
    }
}

size_t sparse_graph::size() const {
    return _N;
}

bool sparse_graph::is_active(uint32_t u) const {
    assert(u < _N);
    return _active.get(u);
}

bool sparse_graph::pi_only(uint32_t u) const {
    assert(u < _N);
    return out_degree_non_pi(u) == 0 && in_degree_non_pi(u) == 0;
}

bool sparse_graph::has_edge(uint32_t u, uint32_t v) const {
    assert(is_active(u) && is_active(v));
    return std::binary_search(std::begin(_out_edges[u]), std::end(_out_edges[u]), v);
}

uint32_t sparse_graph::original_label(uint32_t u) const {
    assert(u < _N);
    return _original_labels[u];
}

const bitvector &sparse_graph::active_vertices() const {
    return _active;
}

void sparse_graph::reactivate_vertex(uint32_t u) {
    assert(!is_active(u));

    for (auto v : in_non_pi(u)) {
        _out_edges_non_pi[v].insert(std::lower_bound(std::begin(_out_edges_non_pi[v]), std::end(_out_edges_non_pi[v]), u), u);
        _out_edges[v].insert(std::lower_bound(std::begin(_out_edges[v]), std::end(_out_edges[v]), u), u);
    }

    for (auto v : out_non_pi(u)) {
        _in_edges_non_pi[v].insert(std::lower_bound(std::begin(_in_edges_non_pi[v]), std::end(_in_edges_non_pi[v]), u), u);
        _in_edges[v].insert(std::lower_bound(std::begin(_in_edges[v]), std::end(_in_edges[v]), u), u);
    }

    for (auto v : pi(u)) {
        if (v == u)
            continue;
        _pi_edges[v].insert(std::lower_bound(std::begin(_pi_edges[v]), std::end(_pi_edges[v]), u), u);
        _out_edges[v].insert(std::lower_bound(std::begin(_out_edges[v]), std::end(_out_edges[v]), u), u);
        _in_edges[v].insert(std::lower_bound(std::begin(_in_edges[v]), std::end(_in_edges[v]), u), u);
    }

    _active.set(u);
}

void sparse_graph::remove_vertex(uint32_t u) {
    assert(is_active(u));
    for (auto v : in_non_pi(u)) {
        _out_edges_non_pi[v].erase(std::lower_bound(std::begin(_out_edges_non_pi[v]), std::end(_out_edges_non_pi[v]), u));
        _out_edges[v].erase(std::lower_bound(std::begin(_out_edges[v]), std::end(_out_edges[v]), u));
    }

    for (auto v : out_non_pi(u)) {
        _in_edges_non_pi[v].erase(std::lower_bound(std::begin(_in_edges_non_pi[v]), std::end(_in_edges_non_pi[v]), u));
        _in_edges[v].erase(std::lower_bound(std::begin(_in_edges[v]), std::end(_in_edges[v]), u));
    }

    for (auto v : pi(u)) {
        if (v == u)
            continue;
        _pi_edges[v].erase(std::lower_bound(std::begin(_pi_edges[v]), std::end(_pi_edges[v]), u));
        _out_edges[v].erase(std::lower_bound(std::begin(_out_edges[v]), std::end(_out_edges[v]), u));
        _in_edges[v].erase(std::lower_bound(std::begin(_in_edges[v]), std::end(_in_edges[v]), u));
    }

    _active.reset(u);
}

void sparse_graph::remove_edge(uint32_t u, uint32_t v) {
    assert(is_active(u) && is_active(v) && has_edge(u, v));
    if (u == v) {
        _pi_edges[u].erase(std::lower_bound(std::begin(_pi_edges[u]), std::end(_pi_edges[u]), u));
    } else if (has_edge(v, u)) {
        _pi_edges[u].erase(std::lower_bound(std::begin(_pi_edges[u]), std::end(_pi_edges[u]), v));
        _pi_edges[v].erase(std::lower_bound(std::begin(_pi_edges[v]), std::end(_pi_edges[v]), u));

        _out_edges_non_pi[v].insert(std::lower_bound(std::begin(_out_edges_non_pi[v]), std::end(_out_edges_non_pi[v]), u), u);
        _in_edges_non_pi[u].insert(std::lower_bound(std::begin(_in_edges_non_pi[u]), std::end(_in_edges_non_pi[u]), v), v);
    } else {
        _out_edges_non_pi[u].erase(std::lower_bound(std::begin(_out_edges_non_pi[u]), std::end(_out_edges_non_pi[u]), v));
        _in_edges_non_pi[v].erase(std::lower_bound(std::begin(_in_edges_non_pi[v]), std::end(_in_edges_non_pi[v]), u));
    }
    _out_edges[u].erase(std::lower_bound(std::begin(_out_edges[u]), std::end(_out_edges[u]), v));
    _in_edges[v].erase(std::lower_bound(std::begin(_in_edges[v]), std::end(_in_edges[v]), u));
}

void sparse_graph::add_edge(uint32_t u, uint32_t v) {
    assert(is_active(u) && is_active(v));
    if (has_edge(u, v))
        return;

    if (u == v) {
        _pi_edges[u].insert(std::lower_bound(std::begin(_pi_edges[u]), std::end(_pi_edges[u]), u), u);
    } else if (has_edge(v, u)) {
        _in_edges_non_pi[u].erase(std::lower_bound(std::begin(_in_edges_non_pi[u]), std::end(_in_edges_non_pi[u]), v));
        _out_edges_non_pi[v].erase(std::lower_bound(std::begin(_out_edges_non_pi[v]), std::end(_out_edges_non_pi[v]), u));

        _pi_edges[u].insert(std::lower_bound(std::begin(_pi_edges[u]), std::end(_pi_edges[u]), v), v);
        _pi_edges[v].insert(std::lower_bound(std::begin(_pi_edges[v]), std::end(_pi_edges[v]), u), u);
    } else {
        _out_edges_non_pi[u].insert(std::lower_bound(std::begin(_out_edges_non_pi[u]), std::end(_out_edges_non_pi[u]), v), v);
        _in_edges_non_pi[v].insert(std::lower_bound(std::begin(_in_edges_non_pi[v]), std::end(_in_edges_non_pi[v]), u), u);
    }
    _out_edges[u].insert(std::lower_bound(std::begin(_out_edges[u]), std::end(_out_edges[u]), v), v);
    _in_edges[v].insert(std::lower_bound(std::begin(_in_edges[v]), std::end(_in_edges[v]), u), u);
}

size_t sparse_graph::out_degree(uint32_t u) const {
    assert(u < _N);
    return _out_edges[u].size();
}

size_t sparse_graph::out_degree_non_pi(uint32_t u) const {
    assert(u < _N);
    return _out_edges_non_pi[u].size();
}

size_t sparse_graph::in_degree(uint32_t u) const {
    assert(u < _N);
    return _in_edges[u].size();
}

size_t sparse_graph::in_degree_non_pi(uint32_t u) const {
    assert(u < _N);
    return _in_edges_non_pi[u].size();
}

size_t sparse_graph::pi_degree(uint32_t u) const {
    assert(u < _N);
    return _pi_edges[u].size();
}

const std::vector<uint32_t> &sparse_graph::out(uint32_t u) const {
    assert(u < _N);
    return _out_edges[u];
}

const std::vector<uint32_t> &sparse_graph::out_non_pi(uint32_t u) const {
    assert(u < _N);
    return _out_edges_non_pi[u];
}

const std::vector<uint32_t> &sparse_graph::in(uint32_t u) const {
    assert(u < _N);
    return _in_edges[u];
}

const std::vector<uint32_t> &sparse_graph::in_non_pi(uint32_t u) const {
    assert(u < _N);
    return _in_edges_non_pi[u];
}

const std::vector<uint32_t> &sparse_graph::pi(uint32_t u) const {
    assert(u < _N);
    return _pi_edges[u];
}

std::istream &operator>>(std::istream &is, sparse_graph &g) {
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

    size_t N, E, t, u, v;
    fscan(N);
    fscan(E);
    fscan(t);

    g._N = N;
    g._active = bitvector(N);
    g._active.fill();
    g._out_edges.resize(N);
    g._out_edges_non_pi.resize(N);
    g._in_edges.resize(N);
    g._in_edges_non_pi.resize(N);
    g._pi_edges.resize(N);
    g._original_labels.resize(N);

    std::iota(std::begin(g._original_labels), std::end(g._original_labels), 0);

    u = 0;
    while (u < N) {
        fscan(v);
        if (v > 0) {
            g._out_edges[u].push_back(v - 1);
            g._in_edges[v - 1].push_back(u);
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
        std::sort(std::begin(g._in_edges[i]), std::end(g._in_edges[i]));
        std::sort(std::begin(g._out_edges[i]), std::end(g._out_edges[i]));
        std::set_intersection(
            std::begin(g._in_edges[i]), std::end(g._in_edges[i]),
            std::begin(g._out_edges[i]), std::end(g._out_edges[i]),
            std::back_inserter(g._pi_edges[i]));
        std::set_difference(
            std::begin(g._in_edges[i]), std::end(g._in_edges[i]),
            std::begin(g._out_edges[i]), std::end(g._out_edges[i]),
            std::back_inserter(g._in_edges_non_pi[i]));
        std::set_difference(
            std::begin(g._out_edges[i]), std::end(g._out_edges[i]),
            std::begin(g._in_edges[i]), std::end(g._in_edges[i]),
            std::back_inserter(g._out_edges_non_pi[i]));
    }
    return is;
}

bool includes(const std::vector<uint32_t> &a, const std::vector<uint32_t> &b) {
    return std::includes(std::begin(a), std::end(a), std::begin(b), std::end(b));
}

uint32_t intersection_size(const std::vector<uint32_t> &a, const std::vector<uint32_t> &b) {
    auto first1 = std::begin(a), last1 = std::end(a);
    auto first2 = std::begin(b), last2 = std::end(b);
    uint32_t count = 0;
    while (first1 != last1 && first2 != last2) {
        if (*first1 < *first2)
            ++first1;
        else if (*first2 < *first1)
            ++first2;
        else {
            ++count;
            ++first1;
            ++first2;
        }
    }
    return count;
}
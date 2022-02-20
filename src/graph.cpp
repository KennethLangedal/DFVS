#include "graph.hpp"

void graph::deactive_edge(size_t u, size_t v) {
    assert(u < N * 64 && v < N * 64 && test(out(u), v) && test(in(v), u));
    reset(out_edges[u], v);
    reset(in_edges[v], u);
}

void graph::deactive_vertex(size_t u) {
    assert(u < N * 64);
    reset(active, u);
    visit(out_edges[u], [&](size_t v) { reset(in_edges[v], u); });
    visit(in_edges[u], [&](size_t v) { reset(out_edges[v], u); });
}

void graph::deactive_single_out(size_t u) {
    assert(u < N * 64 && popcount(out_edges[u]) == 1);
    reset(active, u);
    size_t v;
    visit(out_edges[u], [&](size_t w) {
        reset(in_edges[w], u);
        v = w;
    });
    visit(in_edges[u], [&](size_t w) {
        reset(out_edges[w], u);
        set(out_edges[w], v);
    });
    in_edges[v] |= in_edges[u];
}

void graph::deactive_single_in(size_t u) {
    assert(u < N * 64 && popcount(in_edges[u]) == 1);
    reset(active, u);
    size_t v;
    visit(in_edges[u], [&](size_t w) {
        reset(out_edges[w], u);
        v = w;
    });
    visit(out_edges[u], [&](size_t w) {
        reset(in_edges[w], u);
        set(in_edges[w], v);
    });
    out_edges[v] |= out_edges[u];
}

size_t graph::degree(size_t u) const {
    assert(u < N * 64);
    return popcount(out_edges[u]);
}

size_t graph::in_degree(size_t u) const {
    assert(u < N * 64);
    return popcount(in_edges[u]);
}

bool graph::self_loop(size_t u) const {
    assert(u < N * 64);
    return test(out_edges[u], u);
}

const bitvector<N>& graph::active_vertices() const { return active; }

const bitvector<N>& graph::out(size_t u) const {
    assert(u < N * 64);
    return out_edges[u];
}

const bitvector<N>& graph::in(size_t u) const {
    assert(u < N * 64);
    return in_edges[u];
}

void graph::print_edgelist(std::ostream& os) const {
    visit(active, [&](size_t u) {
        visit(out_edges[u], [&](size_t v) { os << u << " " << v << std::endl; });
    });
}

std::istream& operator>>(std::istream& is, graph& g) {
    size_t N_BUFFER = 1048576;
    std::vector<char> buffer(N_BUFFER);
    size_t buffer_i = N_BUFFER;

    auto get = [&](char& c) {
        if (buffer_i == N_BUFFER) {
            is.read(buffer.data(), N_BUFFER);
            buffer_i = 0;
        }
        c = buffer[buffer_i++];
    };

    char c;
    auto fscan = [&](auto& v) {
        get(c);
        v = 0;
        while (c >= '0' && c <= '9') {
            v = (v * 10) + (c - '0');
            get(c);
        }
    };

    size_t n, E, tmp, u, v;
    fscan(n);
    fscan(E);
    fscan(tmp);
    assert(n <= N * 64);
    g.active.fill(~0ull);
    g.out_edges.resize(N * 64, {});
    g.in_edges.resize(N * 64, {});

    u = 0;
    while (u < n) {
        fscan(v);
        if (v > 0) {
            set(g.out_edges[u], v - 1);
            set(g.in_edges[v - 1], u);
        }
        if (c == '\n') u++;
    }
    return is;
}
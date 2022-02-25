#include "graph.hpp"

void graph::add_edge(size_t u, size_t v) {
    assert(u < N * 64 && v < N * 64 && test(active, u) && test(active, v));
    if (test(out(u), v) || test(in(v), u))
        return;
    log.push_back({action::add_edge, u, v, 0, 0});
    set(out_edges[u], v);
    set(in_edges[v], u);
}

void graph::undo_add_edge(size_t u, size_t v) {
    assert(u < N * 64 && v < N * 64 && test(active, u) && test(active, v) && test(out(u), v) && test(in(v), u));
    reset(out_edges[u], v);
    reset(in_edges[v], u);
}

void graph::remove_edge(size_t u, size_t v) {
    assert(u < N * 64 && v < N * 64 && test(active, u) && test(active, v) && test(out(u), v) && test(in(v), u));
    log.push_back({action::remove_edge, u, v, 0, 0});
    reset(out_edges[u], v);
    reset(in_edges[v], u);
}

void graph::undo_remove_edge(size_t u, size_t v) {
    assert(u < N * 64 && v < N * 64 && test(active, u) && test(active, v) && !test(out(u), v) && !test(in(v), u));
    set(out_edges[u], v);
    set(in_edges[v], u);
}

void graph::remove_vertex(size_t u) {
    assert(u < N * 64 && test(active, u));
    log.push_back({action::remove_vertex, u, 0, 0, 0});
    reset(active, u);
    visit(out_edges[u], [&](size_t v) { reset(in_edges[v], u); });
    visit(in_edges[u], [&](size_t v) { reset(out_edges[v], u); });
}

void graph::undo_remove_vertex(size_t u) {
    assert(u < N * 64 && !test(active, u));
    set(active, u);
    visit(out_edges[u], [&](size_t v) { set(in_edges[v], u); });
    visit(in_edges[u], [&](size_t v) { set(out_edges[v], u); });
}

void graph::fold_neighborhood(size_t u) {
    assert(u < N * 64 && popcount(out(u)) == 2 && out(u) == in(u));
    size_t v[2];
    visit(out_edges[u], [&, i = 0](size_t w) mutable { v[i++] = w; });
    remove_vertex(v[0]);
    remove_vertex(v[1]);
    log.push_back({action::fold_neighborhood, u, v[0], v[1], 0});
    visit(out_edges[v[0]] | out_edges[v[1]], [&](size_t w) {
        if (w != u && w != v[0] && w != v[1]) {
            set(out_edges[u], w);
            set(in_edges[w], u);
        }
    });
    visit(in_edges[v[0]] | in_edges[v[1]], [&](size_t w) {
        if (w != u && w != v[0] && w != v[1]) {
            set(in_edges[u], w);
            set(out_edges[w], u);
        }
    });
}

void graph::undo_fold_neighborhood(size_t u, size_t v1, size_t v2) {
    assert(u < N * 64);
    visit(out_edges[u], [&](size_t w) {
        if (w != v1 && w != v2) {
            reset(out_edges[u], w);
            reset(in_edges[w], u);
        }
    });
    visit(in_edges[u], [&](size_t w) {
        if (w != v1 && w != v2) {
            reset(in_edges[u], w);
            reset(out_edges[w], u);
        }
    });
}

void graph::fold_two_one(size_t u, size_t v1, size_t v2, size_t w) {
    assert(u < N * 64 && popcount(out(u)) == 3 && out(u) == in(u) && test(out(v1), v2) && test(out(v2), v1));
    remove_vertex(u);
    remove_vertex(w);
    visit(out(w), [&](size_t x) {
        if (x != v1 && x != v2 && x != u) {
            add_edge(v1, x);
            add_edge(v2, x);
        }
    });
    visit(in(w), [&](size_t x) {
        if (x != v1 && x != v2 && x != u) {
            add_edge(x, v1);
            add_edge(x, v2);
        }
    });
    log.push_back({action::fold_two_one, u, v1, v2, w});
}

void graph::undo_fold_two_one(size_t u, size_t v1, size_t v2, size_t w) {
    
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

const bitvector<N> &graph::active_vertices() const { return active; }

const bitvector<N> &graph::out(size_t u) const {
    assert(u < N * 64);
    return out_edges[u];
}

const bitvector<N> &graph::in(size_t u) const {
    assert(u < N * 64);
    return in_edges[u];
}

size_t graph::get_timestamp() const {
    return log.size();
}

void graph::unfold_graph(size_t time, bitvector<N> &fvs) {
    while (log.size() > time) {
        auto [t, u, v1, v2, w] = log.back();
        log.pop_back();
        switch (t) {
        case action::add_edge:
            undo_add_edge(u, v1);
            break;
        case action::remove_edge:
            undo_remove_edge(u, v1);
            break;
        case action::remove_vertex:
            undo_remove_vertex(u);
            break;
        case action::fold_neighborhood:
            if (test(fvs, u)) {
                reset(fvs, u);
                set(fvs, v1);
                set(fvs, v2);
            } else {
                set(fvs, u);
            }
            undo_fold_neighborhood(u, v1, v2);
            break;
        case action::fold_two_one:
            if (test(fvs, v1) && test(fvs, v2)) {
                set(fvs, w);
            } else {
                set(fvs, u);
            }
            undo_fold_two_one(u, v1, v2, w);
            break;

        default:
            break;
        }
    }
}

void graph::print_edgelist(std::ostream &os) const {
    visit(active, [&](size_t u) {
        visit(out_edges[u], [&](size_t v) { os << u << " " << v << std::endl; });
    });
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

    char c;
    auto fscan = [&](auto &v) {
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
    if (n > N * 64) {
        std::cout << "Too Large" << std::endl;
        exit(0);
    }
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
        if (c == '\n')
            u++;
    }
    return is;
}
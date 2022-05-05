#include "sparse_local_search.hpp"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <math.h>

void sparse_local_search::insert(uint32_t u, uint32_t i) {
    assert(!config.get(u));
    if (first == N || last == N) { // First element
        first = u;
        last = u;
        next[u] = N;
        prev[u] = N;
    } else if (i == N) { // new end
        next[last] = u;
        prev[u] = last;
        next[u] = N;
        last = u;
    } else {
        if (prev[i] == N) { // new start
            first = u;
            prev[u] = N;
        } else {
            prev[u] = prev[i];
            next[prev[u]] = u;
        }
        prev[i] = u;
        next[u] = i;
    }
    config.set(u);
    assign_label(u);
}

void sparse_local_search::remove(uint32_t u) {
    assert(config.get(u));
    if (prev[u] == N && next[u] == N) { // last element
        first = N;
        last = N;
    } else if (prev[u] == N) { // new start
        first = next[u];
        prev[next[u]] = N;
    } else if (next[u] == N) { // new end
        last = prev[u];
        next[prev[u]] = N;
    } else {
        next[prev[u]] = next[u];
        prev[next[u]] = prev[u];
    }
    config.reset(u);
}

void sparse_local_search::assign_label(uint32_t u) {
    uint64_t low = (u == first) ? std::numeric_limits<uint64_t>::min() : config_order[prev[u]];
    uint64_t high = (u == last) ? std::numeric_limits<uint64_t>::max() : config_order[next[u]];

    if (high == low || high == low + 1) {
        relable();
        assign_label(u);
    } else {
        config_order[u] = low + (high - low) / 2;
    }
}

void sparse_local_search::relable() {
    uint64_t increment = std::numeric_limits<uint64_t>::max() / (config.popcount() + 1);
    uint32_t i = first;
    uint64_t v = std::numeric_limits<uint64_t>::min() + increment;
    while (i != N) {
        config_order[i] = v;
        v += increment;
        i = next[i];
    }
}

void sparse_local_search::search_step(const sparse_graph &g) {
    uint32_t u = dist_int(reng);
    while (config.get(u))
        u = dist_int(reng);
    bool pos = dist_real(reng) > 0.5;
    auto [i, score] = move_score(g, u, pos);
    if (score <= 0 || exp(-score / T) > dist_real(reng)) {
        apply_move(g, u, i, pos);
        if (N - config.popcount() < best_dfvs.popcount()) {
            best_dfvs.set_not(config);
        }
    }
}

std::pair<uint32_t, float> sparse_local_search::move_score(const sparse_graph &g, uint32_t u, bool pos) const {
    float res = -1;
    uint32_t i;
    if (pos) {
        i = N;
        for (auto v : g.out(u)) {
            if (config.get(v) && (i == N || config_order[v] < config_order[i])) {
                i = v;
            }
        }
        if (i != N) {
            for (auto v : g.in(u)) {
                if (config.get(v) && config_order[v] >= config_order[i]) {
                    ++res;
                }
            }
        }
    } else {
        i = first;
        for (auto v : g.in(u)) {
            if (config.get(v) && (i == N || config_order[v] > config_order[i])) {
                i = v;
            }
        }
        if (i != N) {
            for (auto v : g.out(u)) {
                if (config.get(v) && config_order[v] <= config_order[i]) {
                    ++res;
                }
            }
            i = next[i];
        }
    }
    return {i, res};
}

void sparse_local_search::apply_move(const sparse_graph &g, uint32_t u, uint32_t i, bool pos) {
    insert(u, i);
    if (pos) {
        for (auto v : g.in(u)) {
            if (config.get(v) && config_order[v] > config_order[u]) {
                remove(v);
            }
        }
    } else {
        for (auto v : g.out(u)) {
            if (config.get(v) && config_order[v] < config_order[u]) {
                remove(v);
            }
        }
    }
}

sparse_local_search::sparse_local_search(const sparse_graph &g, double T, size_t seed)
    : next(g.size()), prev(g.size()), vertices(g.size()), first(g.size()), last(g.size()), N(g.size()),
      config_order(g.size()), config(g.size()), best_dfvs(g.size()),
      reng(seed), dist_int(0, g.size() - 1), dist_real(0.0, 1.0), T(T) {
    best_dfvs.fill();
    std::iota(std::begin(vertices), std::end(vertices), (uint32_t)0);
}

void sparse_local_search::search(const sparse_graph &g, size_t iterations) {
    for (size_t i = 0; i < iterations; ++i) {
        search_step(g);
    }
}

void sparse_local_search::set_temperature(double T) {
    this->T = T;
}

void sparse_local_search::set_solution(const sparse_graph &g, const bitvector &fvs) {
    config.set_not(fvs);
    shuffle_solution(g);
    if (N - config.popcount() < best_dfvs.popcount()) {
        best_dfvs.set_not(config);
    }
}

void sparse_local_search::shuffle_solution(const sparse_graph &g) {
    std::shuffle(std::begin(vertices), std::end(vertices), reng);
    std::fill(std::begin(config_order), std::end(config_order), 0);

    std::function<void(uint32_t)> visit = [&](uint32_t u) {
        if (!config.get(u) || config_order[u] != 0)
            return;

        config_order[u] = 1;
        for (auto v : g.out(u))
            visit(v);

        if (first == N) {
            last = u;
            first = u;
            next[u] = N;
            prev[u] = N;
        } else {
            prev[first] = u;
            prev[u] = N;
            next[u] = first;
            first = u;
        }
    };

    first = N;

    for (auto u : vertices) {
        visit(u);
    }

    relable();
}

void sparse_local_search::greedy_one_zero_swaps(const sparse_graph &g) {
    bool found = true;
    std::shuffle(std::begin(vertices), std::end(vertices), reng);
    while (found) {
        found = false;
        for (auto u : vertices) {
            if (config.get(u))
                continue;
            auto [i_t, score_t] = move_score(g, u, true);
            if (score_t < 0) {
                apply_move(g, u, i_t, true);
                found = true;
                continue;
            }
            auto [i_f, score_f] = move_score(g, u, false);
            if (score_f < 0) {
                apply_move(g, u, i_f, false);
                found = true;
                continue;
            }
        }
    }
    if (N - config.popcount() < best_dfvs.popcount()) {
        best_dfvs.set_not(config);
    }
}

void sparse_local_search::greedy_two_one_swaps(const sparse_graph &g) {
}

const bitvector &sparse_local_search::get_best() const {
    return best_dfvs;
}

size_t sparse_local_search::get_current() const {
    return N - config.popcount();
}
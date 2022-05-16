#include "local_search_edgew.hpp"
#include <algorithm>
#include <cassert>

void local_search_edgew::_insert(uint32_t u, uint32_t i) {
    assert(_not_config.get(u));
    if (_first == _N || _last == _N) { // First element
        _first = u;
        _last = u;
        _next[u] = _N;
        _prev[u] = _N;
    } else if (i == _N) { // new end
        _next[_last] = u;
        _prev[u] = _last;
        _next[u] = _N;
        _last = u;
    } else {
        if (_prev[i] == _N) { // new start
            _first = u;
            _prev[u] = _N;
        } else {
            _prev[u] = _prev[i];
            _next[_prev[u]] = u;
        }
        _prev[i] = u;
        _next[u] = i;
    }
    _config.set(u);
    _not_config.reset(u);
    _assign_label(u);
}

void local_search_edgew::_remove(uint32_t u) {
    assert(_config.get(u));
    if (_prev[u] == _N && _next[u] == _N) { // last element
        _first = _N;
        _last = _N;
    } else if (_prev[u] == _N) { // new start
        _first = _next[u];
        _prev[_next[u]] = _N;
    } else if (_next[u] == _N) { // new end
        _last = _prev[u];
        _next[_prev[u]] = _N;
    } else {
        _next[_prev[u]] = _next[u];
        _prev[_next[u]] = _prev[u];
    }
    _config.reset(u);
    _not_config.set(u);
}

void local_search_edgew::_assign_label(uint32_t u) {
    uint64_t low = (u == _first) ? std::numeric_limits<uint64_t>::min() : _config_order[_prev[u]];
    uint64_t high = (u == _last) ? std::numeric_limits<uint64_t>::max() : _config_order[_next[u]];

    if (high == low || high == low + 1) {
        _relable();
        _assign_label(u);
    } else {
        _config_order[u] = low + (high - low) / 2;
    }
}

void local_search_edgew::_relable() {
    uint64_t increment = std::numeric_limits<uint64_t>::max() / (_config.popcount() + 1);
    uint32_t i = _first;
    uint64_t v = std::numeric_limits<uint64_t>::min() + increment;
    while (i != _N) {
        _config_order[i] = v;
        v += increment;
        i = _next[i];
    }
}

std::tuple<uint32_t, int32_t> local_search_edgew::_move_score_first_out(uint32_t u) {
    int32_t cost = -1;
    uint32_t i = _N;
    // std::sort(std::begin(_out_edges[u]), std::end(_out_edges[u]), [&](auto a, auto b) {
    //     if (!_config.get(a.first))
    //         return false;
    //     if (!_config.get(b.first))
    //         return true;
    //     return _config_order[a.first] < _config_order[b.first];
    // });

    for (auto [v, e] : _out_edges[u]) {
        if (_config.get(v) && (i == _N || _config_order[v] < _config_order[i])) {
            i = v;
        }
    }
    if (i != _N) {
        for (auto [v, e] : _in_edges[u]) {
            if (_config.get(v) && _config_order[v] >= _config_order[i]) {
                cost += _edge_weights[e];
            }
        }
    }
    return {i, cost};
}

std::tuple<uint32_t, int32_t> local_search_edgew::_move_score_last_in(uint32_t u) {
    int32_t cost = -1;
    uint32_t i = _first;
    // std::sort(std::begin(_in_edges[u]), std::end(_in_edges[u]), [&](auto a, auto b) {
    //     if (!_config.get(a.first))
    //         return false;
    //     if (!_config.get(b.first))
    //         return true;
    //     return _config_order[a.first] > _config_order[b.first];
    // });

    for (auto [v, e] : _in_edges[u]) {
        if (_config.get(v) && (i == _N || _config_order[v] > _config_order[i])) {
            i = v;
        }
    }
    if (i != _N) {
        for (auto [v, e] : _out_edges[u]) {
            if (_config.get(v) && _config_order[v] <= _config_order[i]) {
                cost += _edge_weights[e];
            }
        }
        i = _next[i];
    }
    return {i, cost};
}

void local_search_edgew::_apply_move(uint32_t u, uint32_t i, bool pos) {
    _insert(u, i);
    uint32_t t = 1;
    if (pos) {
        std::shuffle(std::begin(_in_edges[u]), std::end(_in_edges[u]), _reng);
        for (auto [v, e] : _in_edges[u]) {
            if (_config.get(v) && _config_order[v] > _config_order[u]) {
                _edge_weights[e] += t++;
                _remove(v);
            }
        }
    } else {
        std::shuffle(std::begin(_out_edges[u]), std::end(_out_edges[u]), _reng);
        for (auto [v, e] : _out_edges[u]) {
            if (_config.get(v) && _config_order[v] < _config_order[u]) {
                _edge_weights[e] += t++;
                _remove(v);
            }
        }
    }
}

local_search_edgew::local_search_edgew(const sparse_graph &g)
    : _out_edges(g.size()), _in_edges(g.size()), _edge_weights{}, _vertices(g.size()),
      _next(g.size()), _prev(g.size()), _first(g.size()), _last(g.size()), _N(g.size()),
      _config_order(g.size()), _config(g.size()), _not_config(g.size()), _best_dfvs(g.size()),
      _reng(0) {

    _best_dfvs.fill();
    _not_config.fill();

    std::iota(std::begin(_vertices), std::end(_vertices), 0);

    uint32_t edge_counter = 0;
    for (auto u : g.active_vertices()) {
        for (auto v : g.out(u)) {
            _out_edges[u].push_back({v, edge_counter});
            _in_edges[v].push_back({u, edge_counter});
            edge_counter++;
        }
    }
    _edge_weights.resize(edge_counter, 1);
}

void local_search_edgew::search(size_t iterations) {
    for (size_t i = 0; i < iterations; i++) {
        uint32_t best_u, best_i;
        int32_t best_score = std::numeric_limits<int32_t>::max();
        bool best_pos;
        for (auto u : _vertices) {
            if (_config.get(u))
                continue;
            auto [t_i, t_score] = _move_score_first_out(u);
            auto [f_i, f_score] = _move_score_last_in(u);
            if (t_score < best_score) {
                best_score = t_score;
                best_i = t_i;
                best_u = u;
                best_pos = true;
            }
            if (f_score < best_score) {
                best_score = f_score;
                best_i = f_i;
                best_u = u;
                best_pos = false;
            }
        }
        if (best_score != std::numeric_limits<int32_t>::max()) {
            _apply_move(best_u, best_i, best_pos);
            if (_not_config.popcount() < _best_dfvs.popcount()) {
                _best_dfvs.set_not(_config);
            }
        }
    }
}

const bitvector &local_search_edgew::get_best() const {
    return _best_dfvs;
}

const bitvector &local_search_edgew::get_current() const {
    return _not_config;
}

void local_search_edgew::shuffle_solution() {
    std::shuffle(std::begin(_vertices), std::end(_vertices), _reng);
    std::fill(std::begin(_config_order), std::end(_config_order), 0);

    std::function<void(uint32_t)> visit = [&](uint32_t u) {
        if (!_config.get(u) || _config_order[u] != 0)
            return;

        _config_order[u] = 1;
        for (auto [v, e] : _out_edges[u])
            visit(v);

        if (_first == _N) {
            _last = u;
            _first = u;
            _next[u] = _N;
            _prev[u] = _N;
        } else {
            _prev[_first] = u;
            _prev[u] = _N;
            _next[u] = _first;
            _first = u;
        }
    };

    _first = _N;

    for (auto u : _vertices) {
        visit(u);
    }

    _relable();
}
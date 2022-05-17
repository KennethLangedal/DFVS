#include "local_search.hpp"
#include <algorithm>
#include <cassert>

void local_search::_insert(uint32_t u, uint32_t i) {
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
    _current_cost--;
    _not_config.reset(u);
    _assign_label(u);

    for (auto v : _out_edges_not_config[u]) {
        _in_edges_not_config[v].erase(std::remove(std::begin(_in_edges_not_config[v]), std::end(_in_edges_not_config[v]), u), std::end(_in_edges_not_config[v]));
        _in_edges_in_config[v].push_back(u);
        // _in_edges_not_config[v].erase(std::lower_bound(std::begin(_in_edges_not_config[v]), std::end(_in_edges_not_config[v]), u));
        // _in_edges_in_config[v].insert(std::lower_bound(std::begin(_in_edges_in_config[v]), std::end(_in_edges_in_config[v]), u), u);
    }
    for (auto v : _out_edges_in_config[u]) {
        _in_edges_not_config[v].erase(std::remove(std::begin(_in_edges_not_config[v]), std::end(_in_edges_not_config[v]), u), std::end(_in_edges_not_config[v]));
        _in_edges_in_config[v].push_back(u);
        // _in_edges_not_config[v].erase(std::lower_bound(std::begin(_in_edges_not_config[v]), std::end(_in_edges_not_config[v]), u));
        // _in_edges_in_config[v].insert(std::lower_bound(std::begin(_in_edges_in_config[v]), std::end(_in_edges_in_config[v]), u), u);
    }
    for (auto v : _in_edges_not_config[u]) {
        _out_edges_not_config[v].erase(std::remove(std::begin(_out_edges_not_config[v]), std::end(_out_edges_not_config[v]), u), std::end(_out_edges_not_config[v]));
        _out_edges_in_config[v].push_back(u);
        // _out_edges_not_config[v].erase(std::lower_bound(std::begin(_out_edges_not_config[v]), std::end(_out_edges_not_config[v]), u));
        // _out_edges_in_config[v].insert(std::lower_bound(std::begin(_out_edges_in_config[v]), std::end(_out_edges_in_config[v]), u), u);
    }
    for (auto v : _in_edges_in_config[u]) {
        _out_edges_not_config[v].erase(std::remove(std::begin(_out_edges_not_config[v]), std::end(_out_edges_not_config[v]), u), std::end(_out_edges_not_config[v]));
        _out_edges_in_config[v].push_back(u);
        // _out_edges_not_config[v].erase(std::lower_bound(std::begin(_out_edges_not_config[v]), std::end(_out_edges_not_config[v]), u));
        // _out_edges_in_config[v].insert(std::lower_bound(std::begin(_out_edges_in_config[v]), std::end(_out_edges_in_config[v]), u), u);
    }
}

void local_search::_remove(uint32_t u) {
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
    _current_cost++;
    _not_config.set(u);

    for (auto v : _out_edges_not_config[u]) {
        _in_edges_in_config[v].erase(std::remove(std::begin(_in_edges_in_config[v]), std::end(_in_edges_in_config[v]), u), std::end(_in_edges_in_config[v]));
        _in_edges_not_config[v].push_back(u);
        // _in_edges_not_config[v].insert(std::lower_bound(std::begin(_in_edges_not_config[v]), std::end(_in_edges_not_config[v]), u), u);
        // _in_edges_in_config[v].erase(std::lower_bound(std::begin(_in_edges_in_config[v]), std::end(_in_edges_in_config[v]), u));
    }
    for (auto v : _out_edges_in_config[u]) {
        _in_edges_in_config[v].erase(std::remove(std::begin(_in_edges_in_config[v]), std::end(_in_edges_in_config[v]), u), std::end(_in_edges_in_config[v]));
        _in_edges_not_config[v].push_back(u);
        // _in_edges_not_config[v].insert(std::lower_bound(std::begin(_in_edges_not_config[v]), std::end(_in_edges_not_config[v]), u), u);
        // _in_edges_in_config[v].erase(std::lower_bound(std::begin(_in_edges_in_config[v]), std::end(_in_edges_in_config[v]), u));
    }
    for (auto v : _in_edges_not_config[u]) {
        _out_edges_in_config[v].erase(std::remove(std::begin(_out_edges_in_config[v]), std::end(_out_edges_in_config[v]), u), std::end(_out_edges_in_config[v]));
        _out_edges_not_config[v].push_back(u);
        // _out_edges_not_config[v].insert(std::lower_bound(std::begin(_out_edges_not_config[v]), std::end(_out_edges_not_config[v]), u), u);
        // _out_edges_in_config[v].erase(std::lower_bound(std::begin(_out_edges_in_config[v]), std::end(_out_edges_in_config[v]), u));
    }
    for (auto v : _in_edges_in_config[u]) {
        _out_edges_in_config[v].erase(std::remove(std::begin(_out_edges_in_config[v]), std::end(_out_edges_in_config[v]), u), std::end(_out_edges_in_config[v]));
        _out_edges_not_config[v].push_back(u);
        // _out_edges_not_config[v].insert(std::lower_bound(std::begin(_out_edges_not_config[v]), std::end(_out_edges_not_config[v]), u), u);
        // _out_edges_in_config[v].erase(std::lower_bound(std::begin(_out_edges_in_config[v]), std::end(_out_edges_in_config[v]), u));
    }
}

void local_search::_assign_label(uint32_t u) {
    uint64_t low = (u == _first) ? std::numeric_limits<uint64_t>::min() : _config_order[_prev[u]];
    uint64_t high = (u == _last) ? std::numeric_limits<uint64_t>::max() : _config_order[_next[u]];

    if (high == low || high == low + 1) {
        _relable();
        _assign_label(u);
    } else {
        _config_order[u] = low + (high - low) / 2;
    }
}

void local_search::_relable() {
    uint64_t increment = std::numeric_limits<uint64_t>::max() / ((_N - _current_cost) + 1);
    uint32_t i = _first;
    uint64_t v = std::numeric_limits<uint64_t>::min() + increment;
    while (i != _N) {
        _config_order[i] = v;
        v += increment;
        i = _next[i];
    }
}

void local_search::_search_step(const sparse_graph &g) {
    uint32_t u = _dist_int(_reng); //*_not_config.begin().advance(_dist_int(_reng) % _not_config.popcount());
    while (!_not_config.get(u)) {
        u = _dist_int(_reng);
    }
    if (_dist_real(_reng) > 0.5) {
        auto [i, score] = _move_score_first_out(g, u);
        if (score <= 0 || exp(-score / _T) > _dist_real(_reng)) {
            _apply_move(g, u, i, true);
            if (_current_cost < _best_cost) {
                _best_cost = _current_cost;
                _best_dfvs = _not_config;
            }
        }
    } else {
        auto [i, score] = _move_score_last_in(g, u);
        if (score <= 0 || exp(-score / _T) > _dist_real(_reng)) {
            _apply_move(g, u, i, false);
            if (_current_cost < _best_cost) {
                _best_cost = _current_cost;
                _best_dfvs = _not_config;
            }
        }
    }
}

std::tuple<uint32_t, float> local_search::_move_score_first_out(const sparse_graph &g, uint32_t u) const {
    float cost = -1;
    uint32_t i = _N;
    for (auto v : _out_edges_in_config[u]) {
        if (i == _N || _config_order[v] < _config_order[i]) {
            i = v;
        }
    }
    if (i != _N) {
        for (auto v : _in_edges_in_config[u]) {
            if (_config_order[v] >= _config_order[i]) {
                ++cost;
            }
        }
    }
    return {i, cost};
}

std::tuple<uint32_t, float> local_search::_move_score_last_in(const sparse_graph &g, uint32_t u) const {
    float cost = -1;
    uint32_t i = _first;
    for (auto v : _in_edges_in_config[u]) {
        if (i == _N || _config_order[v] > _config_order[i]) {
            i = v;
        }
    }
    if (i != _N) {
        for (auto v : _out_edges_in_config[u]) {
            if (_config_order[v] <= _config_order[i]) {
                ++cost;
            }
        }
        i = _next[i];
    }
    return {i, cost};
}

void local_search::_apply_move(const sparse_graph &g, uint32_t u, uint32_t i, bool pos) {
    _insert(u, i);
    static std::vector<uint32_t> tmp;
    tmp.clear();
    if (pos) {
        std::copy(std::begin(_in_edges_in_config[u]), std::end(_in_edges_in_config[u]), std::back_inserter(tmp));
        for (auto v : tmp) {
            if (_config_order[v] > _config_order[u]) {
                _remove(v);
            }
        }
    } else {
        std::copy(std::begin(_out_edges_in_config[u]), std::end(_out_edges_in_config[u]), std::back_inserter(tmp));
        for (auto v : tmp) {
            if (_config_order[v] < _config_order[u]) {
                _remove(v);
            }
        }
    }
}

local_search::local_search(const sparse_graph &g, double T, size_t seed)
    : _next(g.size()), _prev(g.size()), _vertices(g.size()), _first(g.size()), _last(g.size()), _N(g.size()), _current_cost(g.size()), _best_cost(g.size()),
      _out_edges_in_config(g.size()), _in_edges_in_config(g.size()), _out_edges_not_config(g.size()), _in_edges_not_config(g.size()),
      _config_order(g.size()), _not_config(g.size()), _best_dfvs(g.size()),
      _reng(seed), _dist_int(0, g.size() - 1), _dist_real(0.0, 1.0), _T(T) {
    _best_dfvs.fill();
    _not_config.fill();

    std::iota(std::begin(_vertices), std::end(_vertices), (uint32_t)0);

    for (auto u : g.active_vertices()) {
        for (auto v : g.out(u)) {
            _out_edges_not_config[u].push_back(v);
        }
        for (auto v : g.in(u)) {
            _in_edges_not_config[u].push_back(v);
        }
    }
}

void local_search::search(const sparse_graph &g, size_t iterations) {
    for (size_t i = 0; i < iterations; ++i) {
        _search_step(g);
    }
}

void local_search::set_temperature(double T) {
    _T = T;
}

void local_search::set_solution(const sparse_graph &g, const bitvector &fvs) {
    _not_config = fvs;
    _current_cost = _not_config.popcount();
    shuffle_solution(g);
    if (_current_cost < _best_cost) {
        _best_cost = _current_cost;
        _best_dfvs = _not_config;
    }
}

void local_search::shuffle_solution(const sparse_graph &g) {
    std::shuffle(std::begin(_vertices), std::end(_vertices), _reng);
    std::fill(std::begin(_config_order), std::end(_config_order), 0);

    std::function<void(uint32_t)> visit = [&](uint32_t u) {
        if (_not_config.get(u) || _config_order[u] != 0)
            return;

        _config_order[u] = 1;
        for (auto v : g.out(u))
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

void local_search::greedy_one_zero_swaps(const sparse_graph &g) {
    bool found = true;
    std::shuffle(std::begin(_vertices), std::end(_vertices), _reng);
    while (found) {
        found = false;
        for (auto u : _vertices) {
            if (!_not_config.get(u))
                continue;
            auto [i_t, score_t] = _move_score_first_out(g, u);
            auto [i_f, score_f] = _move_score_last_in(g, u);
            if (score_t < 0) {
                _apply_move(g, u, i_t, true);
                found = true;
            } else if (score_f < 0) {
                _apply_move(g, u, i_f, false);
                found = true;
            } else if (score_t == 0) {
                _apply_move(g, u, i_t, true);
            } else if (score_f == 0) {
                _apply_move(g, u, i_f, false);
            }
        }
    }
    if (_current_cost < _best_cost) {
        _best_cost = _current_cost;
        _best_dfvs = _not_config;
    }
}

const bitvector &local_search::get_best() const {
    return _best_dfvs;
}

const bitvector &local_search::get_current() const {
    return _not_config;
}

const uint32_t &local_search::get_best_cost() const {
    return _best_cost;
}

const uint32_t &local_search::get_current_cost() const {
    return _current_cost;
}
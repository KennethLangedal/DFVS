#include "local_search.hpp"
#include <algorithm>
#include <cassert>

void local_search::_insert(uint32_t u, uint32_t i) {
    assert(!_config[u]);
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
    _config[u] = true;
    _assign_label(u);

    for (auto v : _out_edges_not_config[u]) {
        _in_edges_not_config[v].erase(std::remove(std::begin(_in_edges_not_config[v]), std::end(_in_edges_not_config[v]), u), std::end(_in_edges_not_config[v]));
        _in_edges_in_config[v].push_back(u);
    }
    for (auto v : _out_edges_in_config[u]) {
        _in_edges_not_config[v].erase(std::remove(std::begin(_in_edges_not_config[v]), std::end(_in_edges_not_config[v]), u), std::end(_in_edges_not_config[v]));
        _in_edges_in_config[v].push_back(u);
    }
    for (auto v : _in_edges_not_config[u]) {
        _out_edges_not_config[v].erase(std::remove(std::begin(_out_edges_not_config[v]), std::end(_out_edges_not_config[v]), u), std::end(_out_edges_not_config[v]));
        _out_edges_in_config[v].push_back(u);
    }
    for (auto v : _in_edges_in_config[u]) {
        _out_edges_not_config[v].erase(std::remove(std::begin(_out_edges_not_config[v]), std::end(_out_edges_not_config[v]), u), std::end(_out_edges_not_config[v]));
        _out_edges_in_config[v].push_back(u);
    }
}

void local_search::_remove(uint32_t u) {
    assert(_config[u]);
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
    _config[u] = false;

    for (auto v : _out_edges_not_config[u]) {
        _in_edges_in_config[v].erase(std::remove(std::begin(_in_edges_in_config[v]), std::end(_in_edges_in_config[v]), u), std::end(_in_edges_in_config[v]));
        _in_edges_not_config[v].push_back(u);
    }
    for (auto v : _out_edges_in_config[u]) {
        _in_edges_in_config[v].erase(std::remove(std::begin(_in_edges_in_config[v]), std::end(_in_edges_in_config[v]), u), std::end(_in_edges_in_config[v]));
        _in_edges_not_config[v].push_back(u);
    }
    for (auto v : _in_edges_not_config[u]) {
        _out_edges_in_config[v].erase(std::remove(std::begin(_out_edges_in_config[v]), std::end(_out_edges_in_config[v]), u), std::end(_out_edges_in_config[v]));
        _out_edges_not_config[v].push_back(u);
    }
    for (auto v : _in_edges_in_config[u]) {
        _out_edges_in_config[v].erase(std::remove(std::begin(_out_edges_in_config[v]), std::end(_out_edges_in_config[v]), u), std::end(_out_edges_in_config[v]));
        _out_edges_not_config[v].push_back(u);
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

void local_search::_search_step() {
    uint32_t u = _dist_int(_reng); //*_not_config.begin().advance(_dist_int(_reng) % _not_config.popcount());
    while (_config[u]) {
        u = _dist_int(_reng);
    }

    double roll = _dist_real(_reng);
    if (roll < 0.50) {
        auto [i, score] = _move_score_first_out(u);
        if (score <= 0 || exp(-score / _T) > _dist_real(_reng)) {
            _apply_move(u, i, true);
        }
    } else {
        auto [i, score] = _move_score_last_in(u);
        if (score <= 0 || exp(-score / _T) > _dist_real(_reng)) {
            _apply_move(u, i, false);
        }
    }
    if (_current_cost < _best_cost) {
        _best_cost = _current_cost;
        std::copy(std::begin(_config), std::end(_config), std::begin(_best_config));
    }
}

std::tuple<uint32_t, float> local_search::_move_score_first_out(uint32_t u) const {
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

std::tuple<uint32_t, float> local_search::_move_score_last_in(uint32_t u) const {
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

void local_search::_apply_move(uint32_t u, uint32_t i, bool pos) {
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

void local_search::_greedy_strictly_improving_move(const sparse_graph &g, uint32_t u) {
    static std::vector<uint32_t> current, next;
    static std::vector<bool> visited(_N, false);

    current.clear();
    std::copy(std::begin(_out_edges_in_config[u]), std::end(_out_edges_in_config[u]), std::back_inserter(current));

    std::fill(std::begin(visited), std::end(visited), false);

    bool loop = false;
    while (!current.empty() && !loop) {
        next.clear();

        for (auto v : current) {

            if (std::find(std::begin(_in_edges_in_config[u]), std::end(_in_edges_in_config[u]), v) != std::end(_in_edges_in_config[u])) { // todo improve (bitvector)
                loop = true;
                break;
            }

            for (auto w : _out_edges_in_config[v]) {
                if (visited[w])
                    continue;
                next.push_back(w);
                visited[w] = true;
            }
        }
        std::swap(current, next);
    }

    if (!loop) {
        _insert(u, _N);
        shuffle_solution(g);
    }
}

void local_search::_one_one_candidates(uint32_t u, std::vector<std::vector<uint32_t>> &cand_list) {
    assert(!_config[u]);

    static std::vector<bool> active(_N, true);
    static std::vector<int32_t> degree_out(_N, 0), degree_in(_N, 0);

    uint32_t v = _first;
    while (v != _N) {
        degree_in[v] = _in_edges_in_config[v].size();
        degree_out[v] = _out_edges_in_config[v].size();
        active[v] = true;
        v = _next[v];
    }
    for (auto v : _out_edges_in_config[u]) {
        degree_in[v]++;
    }
    for (auto v : _in_edges_in_config[u]) {
        degree_out[v]++;
    }

    // remove sources
    v = _first;
    while (v != _N) {
        if (active[v]) {
            if (degree_in[v] == 0) {
                active[v] = false;
                for (auto w : _out_edges_in_config[v]) {
                    degree_in[w]--;
                }
            }
        }
        v = _next[v];
    }

    // remove sinks
    v = _last;
    while (v != _N) {
        if (active[v]) {
            if (degree_out[v] == 0) {
                active[v] = false;
                for (auto w : _in_edges_in_config[v]) {
                    degree_out[w]--;
                }
            }
        }
        v = _prev[v];
    }

    // every vertex is now part of cycle
    int32_t leap_edges = 0;
    for (auto v : _out_edges_in_config[u]) {
        if (active[v])
            leap_edges++;
    }

    if (leap_edges == 0) { // u can be added without creating cycles
        return;
    }

    v = _first;
    while (v != _N) {
        if (active[v]) {
            leap_edges -= degree_in[v];
            if (leap_edges == 0) {
                cand_list[v].push_back(u);
            }
            leap_edges += degree_out[v];
        }
        v = _next[v];
    }
}

bool local_search::_two_one_swap_test(const sparse_graph &g, uint32_t u, uint16_t v, uint32_t w) {
    static std::vector<bool> reachable_from_u(_N, false), reachable_to_u(_N, false);
    static std::vector<uint32_t> current, next;

    if (g.has_edge(u, v) && g.has_edge(v, u))
        return true;

    std::fill(std::begin(reachable_from_u), std::end(reachable_from_u), false);
    std::fill(std::begin(reachable_to_u), std::end(reachable_to_u), false);

    current.clear();
    for (auto x : _out_edges_in_config[u]) {
        if (x != w) {
            current.push_back(x);
            reachable_from_u[x] = true;
        }
    }

    while (!current.empty()) {
        next.clear();
        for (auto x : current) {
            for (auto y : _out_edges_in_config[x]) {
                if (!reachable_from_u[y] && y != w) {
                    next.push_back(y);
                    reachable_from_u[y] = true;
                }
            }
        }
        std::swap(current, next);
    }

    current.clear();
    for (auto x : _in_edges_in_config[u]) {
        if (x != w) {
            current.push_back(x);
            reachable_to_u[x] = true;
        }
    }

    while (!current.empty()) {
        next.clear();
        for (auto x : current) {
            for (auto y : _in_edges_in_config[x]) {
                if (!reachable_to_u[y] && y != w) {
                    next.push_back(y);
                    reachable_to_u[y] = true;
                }
            }
        }
        std::swap(current, next);
    }

    bool out_to_u = false;
    for (auto x : _out_edges_in_config[v]) {
        if (x != w && reachable_to_u[x])
            out_to_u = true;
    }
    bool in_from_u = false;
    for (auto x : _in_edges_in_config[v]) {
        if (x != w && reachable_from_u[x])
            in_from_u = true;
    }

    return (out_to_u && in_from_u) || (out_to_u && g.has_edge(u, v)) || (in_from_u && g.has_edge(v, u));
}

local_search::local_search(const sparse_graph &g, double T, size_t seed)
    : _next(g.size()), _prev(g.size()), _vertices(g.size()), _first(g.size()), _last(g.size()), _N(g.size()), _current_cost(g.size()), _best_cost(g.size()),
      _out_edges_in_config(g.size()), _in_edges_in_config(g.size()), _out_edges_not_config(g.size()), _in_edges_not_config(g.size()),
      _config_order(g.size()), _config(g.size(), false), _best_config(g.size(), false),
      _reng(seed), _dist_int(0, g.size() - 1), _dist_real(0.0, 1.0), _T(T) {

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
        _search_step();
    }
}

bool local_search::check_every_two_one_swap(const sparse_graph &g) {
    static std::vector<std::vector<uint32_t>> one_one_candidates(_N);

    for (uint32_t u = 0; u < _N; ++u) {
        one_one_candidates[u].clear();
    }

    for (uint32_t u = 0; u < _N; ++u) {
        if (_config[u])
            continue;
        _one_one_candidates(u, one_one_candidates);
    }

    for (uint32_t u = 0; u < _N; ++u) {
        for (uint32_t i = 0; i < one_one_candidates[u].size(); ++i) {
            for (uint32_t j = i + 1; j < one_one_candidates[u].size(); ++j) {
                if (!_two_one_swap_test(g, one_one_candidates[u][i], one_one_candidates[u][j], u)) {
                    _remove(u);
                    _insert(one_one_candidates[u][i], _N);
                    _insert(one_one_candidates[u][j], _N);
                    shuffle_solution(g);
                    if (_current_cost < _best_cost) {
                        _best_cost = _current_cost;
                        std::copy(std::begin(_config), std::end(_config), std::begin(_best_config));
                    }
                    return true;
                }
            }
        }
    }
    return false;
}

void local_search::set_temperature(double T) {
    _T = T;
}

void local_search::shuffle_solution(const sparse_graph &g) {
    std::shuffle(std::begin(_vertices), std::end(_vertices), _reng);
    std::fill(std::begin(_config_order), std::end(_config_order), 0);

    std::function<void(uint32_t)> visit = [&](uint32_t u) {
        if (!_config[u] || _config_order[u] != 0)
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

void local_search::greedy_one_zero_swaps() {
    bool found = true;
    std::shuffle(std::begin(_vertices), std::end(_vertices), _reng);
    while (found) {
        found = false;
        for (auto u : _vertices) {
            if (_config[u])
                continue;
            auto [i_t, score_t] = _move_score_first_out(u);
            auto [i_f, score_f] = _move_score_last_in(u);
            if (score_t < 0) {
                _apply_move(u, i_t, true);
                found = true;
            } else if (score_f < 0) {
                _apply_move(u, i_f, false);
                found = true;
            } else if (score_t == 0) {
                _apply_move(u, i_t, true);
            } else if (score_f == 0) {
                _apply_move(u, i_f, false);
            }
        }
    }
    if (_current_cost < _best_cost) {
        _best_cost = _current_cost;
        std::copy(std::begin(_config), std::end(_config), std::begin(_best_config));
    }
}

void local_search::greedy_one_zero_swaps_dfs(const sparse_graph &g) {
    bool found = true;
    std::shuffle(std::begin(_vertices), std::end(_vertices), _reng);

    std::vector<uint32_t> current, next;
    std::vector<bool> visited(_N, false);

    while (found) {
        found = false;
        for (auto u : _vertices) {
            if (_config[u])
                continue;

            current.clear();
            std::copy(std::begin(_out_edges_in_config[u]), std::end(_out_edges_in_config[u]), std::back_inserter(current));

            std::fill(std::begin(visited), std::end(visited), false);

            bool loop = false;
            while (!current.empty() && !loop) {
                next.clear();

                for (auto v : current) {

                    if (g.has_edge(v, u)) {
                        loop = true;
                        break;
                    }

                    for (auto w : _out_edges_in_config[v]) {
                        if (visited[w])
                            continue;
                        next.push_back(w);
                        visited[w] = true;
                    }
                }
                std::swap(current, next);
            }

            if (!loop) {
                _insert(u, _N);
                found = true;
            }
        }
        if (found) {
            shuffle_solution(g);
        }
    }
}

bitvector local_search::get_best() const {
    bitvector fvs(_N);
    for (size_t i = 0; i < _N; i++) {
        if (!_best_config[i])
            fvs.set(i);
    }

    return fvs;
}

float local_search::get_average_move_cost(const sparse_graph &g) const {
    double total = 0.0;
    for (auto u : g.active_vertices()) {
        if (_config[u])
            continue;
        auto [i_t, score_t] = _move_score_first_out(u);
        auto [i_f, score_f] = _move_score_last_in(u);
        total += score_t + score_f;
    }
    return total / (2.0 * _N);
}

uint32_t local_search::get_best_cost() const {
    return _best_cost;
}

uint32_t local_search::get_current_cost() const {
    return _current_cost;
}
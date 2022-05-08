#pragma once
#include "bitvector.hpp"
#include <cassert>
#include <random>

template <class graph>
class local_search {
private:
    std::vector<uint32_t> _next, _prev, _vertices;
    uint32_t _first, _last, _N;

    std::vector<uint64_t> _config_order;

    bitvector _config, _not_config, _best_dfvs;

    std::mt19937 _reng;
    std::uniform_int_distribution<uint32_t> _dist_int;
    std::uniform_real_distribution<double> _dist_real;

    double T;

    void _insert(uint32_t u, uint32_t i);
    void _remove(uint32_t u);
    void _assign_label(uint32_t u);
    void _relable();

    void _search_step(const graph &g);
    std::pair<uint32_t, float> _move_score(const graph &g, uint32_t u, bool pos) const;
    void _apply_move(const graph &g, uint32_t u, uint32_t i, bool pos);

public:
    local_search(const graph &g, double T = 0.4, size_t seed = 0);

    void search(const graph &g, size_t iterations);

    void set_temperature(double T);

    void set_solution(const graph &g, const bitvector &fvs);

    void greedy_one_zero_swaps(const graph &g);

    void shuffle_solution(const graph &g);

    const bitvector &get_best() const;

    const bitvector &get_current() const;
};

template <class graph>
void local_search<graph>::_insert(uint32_t u, uint32_t i) {
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
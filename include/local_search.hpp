#pragma once

#include "graph.hpp"
#include <random>

class local_search {
private:
    std::vector<size_t> _pos_in_config;
    bitvector _config_flag, _best_dfvs;
    mutable bitvector tmp0, tmp1;
    std::mt19937 _reng;
    std::uniform_int_distribution<size_t> _dist_int;
    std::uniform_real_distribution<double> _dist_real;
    size_t _best, _N, _config_c;

    size_t l1, l2;
    double T, a;

    bool imp;

    void _search_step(const graph &g);
    std::pair<size_t, double> _move_score(const graph &g, size_t u, bool pos) const;
    void _apply_move(const graph &g, size_t u, size_t i, bool pos);

public:
    local_search(size_t N, size_t seed);

    void search(const graph &g);

    const bitvector &get_best() const;
};